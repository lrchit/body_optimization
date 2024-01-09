
#include <body_optimization.h>

body_optimization::body_optimization() {

  // 读取参数
  YAML::Node config = YAML::LoadFile("../config/nmpc_config.yaml");

  weightQ = SX::eye(3);

  weightQ(0, 0) = config["WeightQ"]["Kpcom"]["x"].as<double>();
  weightQ(1, 1) = config["WeightQ"]["Kpcom"]["y"].as<double>();
  weightQ(2, 2) = config["WeightQ"]["Kpcom"]["z"].as<double>();

  side_sign_x.push_back(1);
  side_sign_x.push_back(1);
  side_sign_x.push_back(-1);
  side_sign_x.push_back(-1);

  side_sign_y.push_back(-1);
  side_sign_y.push_back(1);
  side_sign_y.push_back(-1);
  side_sign_y.push_back(1);
}

body_optimization::~body_optimization() {}

// rpy转为旋转矩阵
SX body_optimization::rpy2rot(SX rpy) {
  SX Rx = SX::eye(3);
  SX Ry = SX::eye(3);
  SX Rz = SX::eye(3);
  SX R = SX::eye(3);

  Rz(0, 0) = cos(rpy(2));
  Rz(0, 1) = -sin(rpy(2));
  Rz(1, 0) = sin(rpy(2));
  Rz(1, 1) = cos(rpy(2));

  Ry(0, 0) = cos(rpy(1));
  Ry(0, 2) = sin(rpy(1));
  Ry(2, 0) = -sin(rpy(1));
  Ry(2, 2) = cos(rpy(1));

  Rx(1, 1) = cos(rpy(0));
  Rx(1, 2) = -sin(rpy(0));
  Rx(2, 1) = sin(rpy(0));
  Rx(2, 2) = cos(rpy(0));

  R = SX::mtimes({Rz, Ry, Rx});

  return R;
}

//创建求解器
void body_optimization::set_nmpc_solver() {

  // Variables，状态和控制输入为变量
  SX obj = SX::sym("obj", 5);
  // Parameters，含有足端位置,yaw
  SX p = SX::sym("p", 12 + 1);

  SX rotation_mat = rpy2rot(SX::vertcat({obj(Slice(0, 2)), p(12)}));

  SX Foot_default = SX::vertcat({dx, dy, dz});

  // 代价函数
  SX f = SX::sym("cost_fun");
  f = 0;
  for (int i = 0; i < 4; ++i) {
    f += SX::mtimes(
        {(Foot_default - p(Slice(i * 3, i * 3 + 3)) + obj(Slice(2, 5)) +
          SX::mtimes({rotation_mat, SX::vertcat({side_sign_x[i] * length,
                                                 side_sign_y[i] * width, 0})}))
             .T(),
         weightQ,
         (Foot_default - p(Slice(i * 3, i * 3 + 3)) + obj(Slice(2, 5)) +
          SX::mtimes(
              {rotation_mat, SX::vertcat({side_sign_x[i] * length,
                                          side_sign_y[i] * width, 0})}))});
  }

  //构建求解器
  SXDict nlp_prob = {{"x", obj}, {"p", p}, {"f", f}};

  string solver_name = "ipopt";
  Dict nlp_opts;
  nlp_opts["expand"] = true;
  nlp_opts["ipopt.max_iter"] = 100;
  nlp_opts["ipopt.linear_solver"] = "mumps";
  nlp_opts["ipopt.print_level"] = 0;
  nlp_opts["print_time"] = 1;
  nlp_opts["ipopt.acceptable_obj_change_tol"] = 1e-6;
  nlp_opts["ipopt.acceptable_tol"] = 1e-4;
  // nlp_opts["ipopt.warm_start_entire_iterate"] = "yes";
  // nlp_opts["ipopt.tol"] = 1e-4;
  // nlp_opts["ipopt.constr_viol_tol"] = 1e-3;
  // nlp_opts["ipopt.fixed_variable_treatment"] = "relax_bounds";

  solver = nlpsol("nlpsol", solver_name, nlp_prob, nlp_opts);

  // // Generate C code for the NLP functions
  // solver.generate_dependencies("nlp");

  // // Just-in-time compilation?
  // bool jit = false;
  // if (jit) {
  //   // Create a new NLP solver instance using just-in-time compilation
  //   solver = nlpsol("solver", solver_name, "nlp.c");
  // } else {
  //   // Compile the c-code
  //   int flag = system("gcc -fPIC -shared -O3 nlp.c -o nlp.so");
  //   casadi_assert(flag == 0, "Compilation failed");

  //   // Create a new NLP solver instance from the compiled code
  //   solver = nlpsol("solver", solver_name, "nlp.so");
  // }
}

void body_optimization::opti_solution(
    Eigen::Matrix<double, 6, 1> current_states,
    Eigen::Matrix<double, 12, 1> p_feet) {

  //求解参数设置
  // Initial guess and bounds for the optimization variablese
  vector<double> x0;
  x0.push_back(current_states(0));
  x0.push_back(current_states(1));
  for (int i = 0; i < 3; ++i) {
    x0.push_back(current_states(3 + i));
  }

  // Nonlinear bounds
  vector<double> lbx;
  vector<double> ubx;
  double pFoot_x_min =
      fmin(fmin(fmin(fmin(fmin(p_feet(0) + p_feet(3), p_feet(0) + p_feet(6)),
                          p_feet(0) + p_feet(9)),
                     p_feet(3) + p_feet(6)),
                p_feet(3) + p_feet(9)),
           p_feet(6) + p_feet(9));
  double pFoot_x_max =
      fmax(fmax(fmax(fmax(fmax(p_feet(0) + p_feet(3), p_feet(0) + p_feet(6)),
                          p_feet(0) + p_feet(9)),
                     p_feet(3) + p_feet(6)),
                p_feet(3) + p_feet(9)),
           p_feet(6) + p_feet(9));
  double pFoot_y_min =
      fmin(fmin(fmin(fmin(fmin(p_feet(1) + p_feet(4), p_feet(1) + p_feet(7)),
                          p_feet(1) + p_feet(10)),
                     p_feet(4) + p_feet(7)),
                p_feet(4) + p_feet(10)),
           p_feet(7) + p_feet(10));
  double pFoot_y_max =
      fmax(fmax(fmax(fmax(fmax(p_feet(1) + p_feet(4), p_feet(1) + p_feet(7)),
                          p_feet(1) + p_feet(10)),
                     p_feet(4) + p_feet(7)),
                p_feet(4) + p_feet(10)),
           p_feet(7) + p_feet(10));

  lbx.push_back(-0.8);
  lbx.push_back(-0.8);
  lbx.push_back(pFoot_x_min);
  lbx.push_back(pFoot_y_min);
  lbx.push_back(0.2);
  ubx.push_back(0.8);
  ubx.push_back(0.8);
  ubx.push_back(pFoot_x_max);
  ubx.push_back(pFoot_y_max);
  ubx.push_back(0.36);

  // Original parameter values
  vector<double> p0;
  for (int leg = 0; leg < 4; leg++) {
    for (int i = 0; i < 3; i++) {
      p0.push_back(p_feet(i + 3 * leg));
    }
  }
  p0.push_back(current_states(2));

  arg["lbx"] = lbx;
  arg["ubx"] = ubx;
  arg["x0"] = x0;
  arg["p"] = p0;
  res = solver(arg);

  std::cout << "Objective: " << res.at("f") << std::endl;
  // cout << "Optimal solution for p = \n" << arg.at("p") << ":" << endl;
  // cout << "Primal solution: \n" << res.at("x") << endl;
  // cout << "Dual solution (x): " << res.at("lam_x") << endl;
  // cout << "Dual solution (g): " << res.at("lam_g") << endl;

  vector<double> res_all(res.at("x"));
  for (int i = 0; i < 5; i++) {
    opt_obj(i) = res_all.at(i);
  }
}

int main(int argc, char *argv[]) {

  // 足底位置
  Eigen::Matrix<double, 12, 1> p_feet;
  p_feet << 0.183, -0.13205, 0, 0.183, 0.13205, 0, -0.183, -0.13205, 0, -0.183,
      0.13205, 0;
  // std::cout << "p_feet\n" << p_feet << std::endl;

  // 当前状态，参考状态序列，参考控制序列
  Eigen::Matrix<double, 6, 1> current_states;
  current_states << 0, 0, 0, 0, 0, 0.28;

  // 生成nmpc类
  body_optimization nmpc_controller;
  nmpc_controller.set_nmpc_solver();

  // 调用优化
  double iteration_steps = 100;
  for (int i = 0; i < iteration_steps + 1; i++) {

    nmpc_controller.opti_solution(current_states, p_feet);
    //获取最优控制向量
    Eigen::Matrix<double, 5, 1> obj;
    obj = nmpc_controller.get_opt_obj();

    // 更新状态和参考控制序列
    // 取下一次调用mpc时的状态和算出的这次的控制输入
    current_states.segment(0, 2) = obj.segment(0, 2);
    current_states.segment(3, 3) = obj.segment(2, 3);

    std::cout << "\n********* iteration_steps:" << i << " *********"
              << std::endl;
    std::cout << "current_states\n" << current_states.transpose() << std::endl;
  }
}
