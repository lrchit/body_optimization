cmake_minimum_required(VERSION 3.16)
project(body_optimization)
set(CMAKE_CXX_FLAGS "-std=c++11")

link_directories("./lib")
include_directories("./include")
add_executable(body_optimization
  src/body_optimization.cpp
)
target_link_libraries(body_optimization libyaml-cpp.a)
target_link_libraries(body_optimization /usr/local/lib/libcasadi.so.3.7)
