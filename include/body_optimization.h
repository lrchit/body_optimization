
#ifndef BODY_OPTIMIZATION_H
#define BODY_OPTIMIZATION_H

#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>

using namespace casadi;
using namespace std;

class body_optimization {
public:
  body_optimization();
  ~body_optimization();
  //定义求解器
  void set_nmpc_solver();
  //优化求解
  void opti_solution(Eigen::Matrix<double, 6, 1> current_states,
                     Eigen::Matrix<double, 12, 1> p_feet);
  //获取解
  Eigen::Matrix<double, 5, 1> get_opt_obj() { return opt_obj; }

private:
  SX weightQ;

  double dx = 0;
  double dy = 0;
  double dz = -0.28;
  double length = 0.183;
  double width = 0.13205;
  vector<double> side_sign_x;
  vector<double> side_sign_y;

  Eigen::Matrix<double, 5, 1> opt_obj;

  Function solver;     //求解器
  map<string, DM> res; //求解结果
  map<string, DM> arg; //求解参数

  SX rpy2rot(SX rpy);
};

#endif // BODY_OPTIMIZATION_H