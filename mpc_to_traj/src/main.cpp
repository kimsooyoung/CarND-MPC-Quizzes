#include <matplotlibcpp.h>
#include <Eigen/QR>
#include <vector>
#include "traj_helper.h"
#include "helpers.h"
#include "MPC.h"

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

  const double pi = M_PI;

  const int traj_sampling_num = 100;
  const int window_size = 10;
  const int iters = 100;
  const double dt = 0.1;

  // initial trajectory

  // Case 1. line trajectory
  auto traj_points = GetTrajPointsLine(0, window_size, traj_sampling_num);
  // Case 2. circular trajectory
  // auto traj_points = GetTrajPointsCirc(0, window_size, traj_sampling_num);
  
  auto traj_x = std::get<0>(traj_points);
  auto traj_y = std::get<1>(traj_points);

  VectorXd x_veh(window_size);
  VectorXd y_veh(window_size);

  for(long unsigned int i = 0; i < window_size; ++i)
  {
    x_veh[i] = traj_x[i];
    y_veh[i] = traj_y[i];
  }

  // 3rd order polynomial
  auto coeffs = polyfit(x_veh, y_veh, 1);

  // initial state
  double x = 0.0;
  double y = 0.0;
  // double y = 2.0;
  double psi = 0.0;
  double v = 0.0;
  double w = 0.0;
  // We'll use robot frame, therefore cte = 0 is our goal.
  double cte = polyeval(coeffs, x);
  // linear approx for epsi
  double epsi = atan(coeffs[1]);

  // prepare parameters
  map<string, double> mpc_params;
  mpc_params["STEPS"] = window_size;
  // mpc_params["REF_V"] = (2 * pi) / (dt * iters);
  mpc_params["REF_V"] = 1.0;
  mpc_params["DT"] = dt;

  VectorXd state(7);
  state << x, y, psi, v, w, cte, epsi;

  // create mpc instance
  MPC mpc;
  mpc.LoadParams(mpc_params);

  // output container
  vector<double> x_vals = {state[0]};
  vector<double> y_vals = {state[1]};
  vector<double> psi_vals = {state[2]};
  vector<double> v_vals = {state[3]};
  vector<double> w_vals = {state[4]};
  vector<double> cte_vals = {state[5]};
  vector<double> epsi_vals = {state[6]};
  vector<double> a_vals;
  vector<double> alpha_vals;

  for (size_t i = 0; i < 1; ++i) {
    // cout << "Iteration " << i << endl;

    auto vars = mpc.Solve(state, coeffs);

    // cout << "x = " << vars[0] << endl;
    // cout << "y = " << vars[1] << endl;
    // cout << "psi = " << vars[2] << endl;
    // cout << "v = " << vars[3] << endl;
    // cout << "w = " << vars[4] << endl;
    // cout << "cte = " << vars[5] << endl;
    // cout << "epsi = " << vars[6] << endl;
    // cout << "a = " << vars[7] << endl;
    // cout << "alpha = " << vars[8] << endl;
    // cout << endl;

    auto cur_x = vars[0];
    auto cur_y = vars[1];
    auto cur_psi = vars[2];
    auto cur_v = vars[3];
    auto cur_w = vars[4];
    auto cur_cte = vars[5];
    auto cur_epsi = vars[6];
    auto cur_a = vars[7];
    auto cur_alpha = vars[8];

    // x_vals.push_back(cur_x + cur_v * cos(cur_psi) * dt);
    // y_vals.push_back(cur_y + cur_v * sin(cur_psi) * dt);
    // psi_vals.push_back(cur_psi + cur_w * dt);

    x_vals.push_back(cur_x);
    y_vals.push_back(cur_y);
    psi_vals.push_back(cur_psi + cur_w * dt);
    v_vals.push_back(cur_v);
    w_vals.push_back(cur_w);
    cte_vals.push_back(cur_cte);
    epsi_vals.push_back(cur_epsi);
    a_vals.push_back(cur_a);
    alpha_vals.push_back(cur_alpha);

    // update trajectory
    // auto traj_points = GetTrajPointsCirc(i, window_size, traj_sampling_num);
    auto traj_points = GetTrajPointsLine(i, window_size, traj_sampling_num);
    auto traj_x = std::get<0>(traj_points);
    auto traj_y = std::get<1>(traj_points);

    for(long unsigned int i = 0; i < window_size; ++i)
    {
      x_veh[i] = traj_x[i];
      y_veh[i] = traj_y[i];
    }

    // update coeffs
    coeffs = polyfit(x_veh, y_veh, 1);

    const double cte  = polyeval(coeffs, 0.0);
    const double epsi = atan(coeffs[1]);

    // update state
    // state << vars[0], vars[1], vars[2], vars[3], vars[4], cte, epsi;
    auto new_x = cur_x + cur_v * cos(cur_psi) * dt;
    auto new_y = cur_y + cur_v * sin(cur_psi) * dt;
    auto new_psi = cur_psi + cur_w * dt;
    auto new_v = cur_v + cur_a * dt;
    auto new_w = cur_w + cur_alpha * dt;
    
    state << new_x, new_y, new_psi, new_v, new_w, cte, epsi;
  }

  // vector<double> gt_x(traj_sampling_num + 1);
  // vector<double> gt_y(traj_sampling_num + 1);

  // for(int i = 0; i < traj_sampling_num + 1; ++i)
  // {
  //     gt_x[i] = cos(2 * pi * i / traj_sampling_num);
  //     gt_y[i] = sin(2 * pi * i / traj_sampling_num) + 1;
  // }

  auto gt_traj = GetTrajPointsLine(0, traj_sampling_num, traj_sampling_num);
  auto gt_x = std::get<0>(gt_traj);
  auto gt_y = std::get<1>(gt_traj);

  if (false){
    plt::figure(1);
    plt::subplot(3, 1, 1);
    plt::title("X Values");
    plt::plot(x_vals);
    plt::subplot(3, 1, 2);
    plt::title("Y Values");
    plt::plot(y_vals);
    plt::subplot(3, 1, 3);
    plt::title("PSI Values");
    plt::plot(psi_vals);

    plt::figure(2);
    plt::subplot(2, 1, 1);
    plt::title("V");
    plt::plot(v_vals);
    plt::subplot(2, 1, 2);
    plt::title("W");
    plt::plot(w_vals);

    plt::figure(3);
    plt::subplot(2, 1, 1);
    plt::title("CTE");
    plt::plot(cte_vals);
    plt::subplot(2, 1, 2);
    plt::title("EPSI");
    plt::plot(epsi_vals);

    plt::figure(4);
    plt::subplot(2, 1, 1);
    plt::title("Acc");
    plt::plot(a_vals);
    plt::subplot(2, 1, 2);
    plt::title("Anaugular Acc");
    plt::plot(alpha_vals);

    plt::figure(5);
    plt::plot(gt_x, gt_y, "r--"); //plot the x,y
    plt::plot(x_vals, y_vals); //plot the x,y
    plt::grid(true); //show grid

    // plt::subplot(3, 1, 1);
    // plt::title("CTE");
    // plt::plot(cte_vals);
    // plt::subplot(3, 1, 2);
    // plt::title("Delta (Radians)");
    // plt::plot(delta_vals);
    // plt::subplot(3, 1, 3);
    // plt::title("Accel m/s^2");
    // plt::plot(a_vals);

    plt::show();
  }

  return 0;
}