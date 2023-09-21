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
  const double ref_v = 0.2; // v_xy

  // int traj_sampling_num = int(iters * dt);
  // int traj_sampling_num = 100;

  // Retrieve Full Trajectory
  // 10m를 200번만에 간다고 하자.
  auto traj_points = GetTrajPointsSquareShape(0, iters, 40.0, iters);
  auto full_traj_x = std::get<0>(traj_points);
  auto full_traj_y = std::get<1>(traj_points);

  // initial state 선언
  double x = 0.0;
  double y = 0.0;
  double psi = 0.0;
  double v_x = 0.0;
  double v_y = 0.0;
  double w = 0.0;

  // initial trajectory transformation
  vector<double> robot_frame_x_traj(window_size);
  vector<double> robot_frame_y_traj(window_size);

  double cospsi = cos(psi); // 1
  double sinpsi = sin(psi); // 0

  for(long unsigned int j = 0; j < window_size; j++) 
  {
      robot_frame_x_traj[j] = full_traj_x[0 + j] * cospsi + full_traj_y[0 + j] * sinpsi - x * cospsi - y * sinpsi;
      robot_frame_y_traj[j] = -full_traj_x[0 + j] * sinpsi + full_traj_y[0 + j] * cospsi + x * sinpsi - y * cospsi;
  }

  // for(auto x : robot_frame_x_traj)
  //   std::cout << x << " ";
  // std::cout << std::endl;

  // for(auto y : robot_frame_y_traj)
  //   std::cout << y << " ";
  // std::cout << std::endl;

  double cte_x = robot_frame_x_traj[0];
  double cte_y = robot_frame_x_traj[0];
  
  // MPC 구현
  VectorXd state(8);
  state << x, y, psi, v_x, v_y, w, cte_x, cte_y;

  // prepare mpc parameters
  map<string, double> mpc_params;
  mpc_params["STEPS"] = window_size;
  mpc_params["REF_V"] = ref_v;
  mpc_params["DT"] = dt;
  
  mpc_params["MIN_ACC_X"] = -5.0;
  mpc_params["MAX_ACC_X"] = 5.0;
  mpc_params["MIN_ACC_Y"] = -5.0;
  mpc_params["MAX_ACC_Y"] = 5.0;
  mpc_params["MIN_ANG_ACC"] = -3.0;
  mpc_params["MAX_ANG_ACC"] = 3.0;

  mpc_params["W_CTE_X"] = 450.0;
  mpc_params["W_CTE_Y"] = 450.0;
  mpc_params["W_V"] = 20.0;
  mpc_params["W_AX"] = 1.0;
  mpc_params["W_AY"] = 1.0;
  mpc_params["W_ALPHA"] = 1.0;
  mpc_params["W_DELTA_AX"] = 1.0;
  mpc_params["W_DELTA_AY"] = 1.0;
  mpc_params["W_DELTA_ALPHA"] = 1.0;

  // create mpc instance
  MPCQD mpc;
  mpc.LoadParams(mpc_params);

  // parse mpc result
  auto robot_traj_tup = std::make_tuple(robot_frame_x_traj, robot_frame_y_traj);
  
  // output container
  double x_world = 0;
  double y_world = 0;
  double psi_world = 0;
  vector<double> x_vals = {state[0]};
  vector<double> y_vals = {state[1]};

  vector<double> cur_x_traj(window_size);
  vector<double> cur_y_traj(window_size);

  // auto vars = mpc.Solve(state, robot_traj_tup);
  for(long unsigned int i = 0; i < iters; ++i)
  {
    auto vars = mpc.Solve(state, robot_traj_tup);

    auto cur_x = vars[0];
    auto cur_y = vars[1];
    auto cur_psi = vars[2];
    auto cur_vx = vars[3];
    auto cur_vy = vars[4];
    auto cur_w = vars[5];

    x_world += cur_x * cos(psi_world) - cur_y * sin(psi_world);
    y_world += cur_x * sin(psi_world) + cur_y * cos(psi_world);
    psi_world += cur_psi;

    // renew trajectory
    if(i + window_size >= iters){
      int overlap_num = i + window_size - iters + 1;

      for(int j = 0; j < window_size - overlap_num; j++)
      {
        cur_x_traj[j] = full_traj_x[i + j];
        cur_y_traj[j] = full_traj_y[i + j];
      }
      auto last_x_traj = cur_x_traj[window_size - overlap_num - 1];
      auto last_y_traj = cur_y_traj[window_size - overlap_num - 1];

      for(int j = 0; j < overlap_num; j++)
      {
        cur_x_traj[window_size - overlap_num + j] = last_x_traj;
        cur_y_traj[window_size - overlap_num + j] = last_y_traj;
      }
    }
    else{
      for(long unsigned int j = 0; j < window_size; j++)
      {
        cur_x_traj[j] = full_traj_x[i + j];
        cur_y_traj[j] = full_traj_y[i + j];
      }
    }
    
    // Robot Frame conversion cur_x_traj => robot_frame_x_traj
    cospsi = cos(psi_world);
    sinpsi = sin(psi_world);
    
    for(long unsigned int j = 0; j < window_size; j++) 
    {
        robot_frame_x_traj[j] = cur_x_traj[j] * cospsi + cur_y_traj[j] * sinpsi - x_world * cospsi - y_world * sinpsi;
        robot_frame_y_traj[j] = -cur_x_traj[j] * sinpsi + cur_y_traj[j] * cospsi + x_world * sinpsi - y_world * cospsi;
    }

    robot_traj_tup = std::make_tuple(robot_frame_x_traj, robot_frame_y_traj);

    auto cur_cte_x = robot_frame_x_traj[0] - cur_x;
    auto cur_cte_y = robot_frame_y_traj[0] - cur_y;

    std::cout << "cur_x : " << cur_x << std::endl;
    std::cout << "cur_y : " << cur_y << std::endl;
    std::cout << "cur_psi : " << cur_psi << std::endl;
    std::cout << "cur_vx : " << cur_vx << std::endl;
    std::cout << "cur_vy : " << cur_vy << std::endl;
    std::cout << "cur_w : " << cur_w << std::endl;
    std::cout << "cur_cte_x : " << cur_cte_x << std::endl;
    std::cout << "cur_cte_y : " << cur_cte_y << std::endl << std::endl;

    state << 0, 0, 0, cur_vx, cur_vy, cur_w, cur_cte_x, cur_cte_y;

    x_vals.push_back(x_world);
    y_vals.push_back(y_world);

    plt::xlim(-1, 15);
    plt::ylim(-15, 1);
    plt::grid(true); //show grid
    plt::plot(full_traj_x, full_traj_y, "b");
    plt::plot(x_vals, y_vals, "r"); //plot the x,y
    plt::pause(0.05);
  }

  plt::figure(1);
  plt::xlim(-1, 15);
  plt::ylim(-15, 1);
  plt::grid(true); //show grid
  plt::plot(full_traj_x, full_traj_y, "b");
  plt::plot(x_vals, y_vals, "r"); //plot the x,y

  plt::show();
  
  return 0;
}