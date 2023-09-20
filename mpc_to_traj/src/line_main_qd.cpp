#include <matplotlibcpp.h>
#include <Eigen/QR>
#include <vector>
#include "traj_helper.h"
#include "helpers.h"
#include "MPCQD.h"

namespace plt = matplotlibcpp;

using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::string;

// For every steps
// 1. get traj
// 2. convert traj into robot frame
// 3. get coeffs
// 4. solve mpc
// 5. update states
// 6. plot results

int main() {

  // const double pi = M_PI;

  const int window_size = 20;
  const int iters = 200; // 200
  const double dt = 0.1;
  const double ref_v = 1.5; // v_xy

  // int traj_sampling_num = int(iters * dt);
  int traj_sampling_num = 100;

  // Retrieve Full Trajectory
  // 10m를 200번만에 간다고 하자.
  auto traj_points = GetTrajPointsLineWithOption(0, iters, 10.0, iters);
  auto full_traj_x = std::get<0>(traj_points);
  auto full_traj_y = std::get<1>(traj_points);

  vector<double> cur_x_traj(window_size);
  vector<double> cur_y_traj(window_size);  

  for(long unsigned int j = 0; j < window_size; j++)
  {
    cur_x_traj[j] = full_traj_x[0 + j];
    cur_y_traj[j] = full_traj_y[0 + j];
  }

  // initial state 선언
  double x = 0.0;
  double y = 0.0;
  double psi = 0.0;
  double v_x = 0.0;
  double v_y = 0.0;
  double w = 0.0;
  double cte_x = 0.0;
  double cte_y = 0.0;
  
  // trajectory를 robot frame으로 변환
  // initial trajectory transformation
  double cospsi = cos(psi);
  double sinpsi = sin(psi);
  
  vector<double> robot_frame_x_traj(window_size);
  vector<double> robot_frame_y_traj(window_size);

  for(long unsigned int i = 0; i < window_size; ++i)
  {
    const double dx = cur_x_traj[i] - x;
    const double dy = cur_y_traj[i] - y;

    robot_frame_x_traj[i] = dx * cospsi + dy * sinpsi;
    robot_frame_y_traj[i] = dy * cospsi - dx * sinpsi;
  }
  
  // MPC 구현
  VectorXd state(8);
  state << x, y, psi, v_x, v_y, w, cte_x, cte_y;

  // create mpc instance
  MPCQD mpc;
  // TODO Parameter Setup
  // mpc.LoadParams(mpc_params);

  // parse mpc result
  auto cur_traj_tup = std::make_tuple(robot_frame_x_traj, robot_frame_y_traj);
  auto vars = mpc.Solve(state, cur_traj_tup);

  // // TODO: 남은 trajectory 수가 window size보다 작을 때 처리
  // for(long unsigned int i = 0; i < 50; i++){

  //   // renew trajectory
  //   for(long unsigned int j = 0; j < window_size; j++)
  //   {
  //     cur_x_traj[j] = full_traj_x[i+j];
  //     cur_y_traj[j] = full_traj_y[i+j];
  //   }

  //   // parse mpc result
  //   auto cur_traj_tup = std::make_tuple(x, y);
  //   auto vars = mpc.Solve(state, cur_traj_tup);

  //   // update state


  //   plt::clf();
  //   plt::xlim(-1, 50);
  //   plt::grid(true); //show grid
  //   plt::plot(cur_x_traj, cur_y_traj, "r"); //plot the x,y
  //   plt::pause(0.01);

  // }

  // plt::figure(1);
  // plt::show();

  return 0;
}