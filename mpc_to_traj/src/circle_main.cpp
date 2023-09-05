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
using std::sin;
using std::cos;

// For every steps
// 1. get traj
// 2. convert traj into robot frame
// 3. get coeffs
// 4. solve mpc
// 5. update states
// 6. plot results

int main() {

  const double pi = M_PI;

  const int window_size = 9;
  const int iters = 200; // 200
  const double dt = 0.1;

  int traj_sampling_num = int(iters * dt);

  // Case 2. circular trajectory
  auto traj_points = GetTrajPointsCirc(0, window_size, iters);
  auto traj_x = std::get<0>(traj_points);
  auto traj_y = std::get<1>(traj_points);

  // initial state
  double x = 0.0;
  double y = 0.0;
  double psi = 0.0;
  double v = 0.0;
  double w = 0.0;
  // We'll use robot frame, therefore cte = 0 is our goal.

  VectorXd x_veh(window_size);
  VectorXd y_veh(window_size);
  double cospsi = cos(psi);
  double sinpsi = sin(psi);

  // initial trajectory
  for(long unsigned int i = 0; i < window_size; ++i)
  {
    const double dx = traj_x[i] - x;
    const double dy = traj_y[i] - y;
    x_veh[i] = dx * cospsi + dy * sinpsi;
    y_veh[i] = dy * cospsi - dx * sinpsi;
  }

  // 1rd order polynomial
  auto coeffs = polyfit(x_veh, y_veh, 1);
  // 3rd order polynomial
  // auto coeffs = polyfit(x_veh, y_veh, 3);

  // double cte = polyeval(coeffs, x) - y;
  double cte = polyeval(coeffs, 0.0);
  // linear approx for epsi
  // double epsi = psi - atan(coeffs[1]);
  double epsi = atan(coeffs[1]);

  // prepare parameters
  map<string, double> mpc_params;
  mpc_params["STEPS"] = window_size;
  mpc_params["REF_V"] = 1 * (2 * pi) / (dt * iters);
  // mpc_params["REF_V"] = 1;
  mpc_params["DT"] = dt;
  mpc_params["MIN_ACC"] = -5.0;
  mpc_params["MAX_ACC"] = 5.0;
  mpc_params["MIN_ANG_ACC"] = -3.0;
  mpc_params["MAX_ANG_ACC"] = 3.0;

  mpc_params["W_CTE"] = 100.0;
  mpc_params["W_EPSI"] = 100.0;
  mpc_params["W_V"] = 10.0;
  mpc_params["W_A"] = 1.0;
  mpc_params["W_ALPHA"] = 1.0;
  mpc_params["W_DELTA_A"] = 0.0;
  mpc_params["W_DELTA_ALPHA"] = 0.0;

  VectorXd state(7);
  state << x, y, psi, v, w, cte, epsi;

  // create mpc instance
  MPC mpc;
  mpc.LoadParams(mpc_params);

  // output container
  double x_world = 0;
  double y_world = 0;
  double psi_world = 0;
  vector<double> x_vals = {state[0]};
  vector<double> y_vals = {state[1]};
  vector<double> psi_vals = {state[2]};
  vector<double> v_vals = {state[3]};
  vector<double> w_vals = {state[4]};
  vector<double> cte_vals = {state[5]};
  vector<double> epsi_vals = {state[6]};
  vector<double> a_vals;
  vector<double> alpha_vals;

  // auto vars = mpc.Solve(state, coeffs);
  // return 0;

  // for (size_t i = 0; i < iters; ++i) {
  for (size_t i = 0; i < 100; ++i) {

    auto vars = mpc.Solve(state, coeffs);

    auto cur_x = vars[0];
    auto cur_y = vars[1];
    auto cur_psi = vars[2];
    auto cur_v = vars[3];
    auto cur_w = vars[4];
    auto cur_cte = vars[5];
    auto cur_epsi = vars[6];
    auto cur_a = vars[7];
    auto cur_alpha = vars[8];

    x_world += cur_x * cos(psi_world);
    y_world += cur_x * sin(psi_world);
    psi_world += cur_psi;

    cout << "x = " << x_world << endl;
    cout << "y = " << y_world << endl;
    cout << "psi = " << psi_world << endl;
    cout << "v = " << vars[3] << endl;
    cout << "w = " << vars[4] << endl;
    cout << "cte = " << vars[5] << endl;
    cout << "epsi = " << vars[6] << endl;
    cout << "a = " << vars[7] << endl;
    cout << "alpha = " << vars[8] << endl;
    cout << endl;

    x_vals.push_back(x_world);
    y_vals.push_back(y_world);
    psi_vals.push_back(psi_world);
    v_vals.push_back(cur_v);
    w_vals.push_back(cur_w);
    cte_vals.push_back(cur_cte);
    epsi_vals.push_back(cur_epsi);
    a_vals.push_back(cur_a);
    alpha_vals.push_back(cur_alpha);

    // get new trajectory
    traj_points = GetTrajPointsCirc(i, window_size, iters);
    traj_x = std::get<0>(traj_points);
    traj_y = std::get<1>(traj_points);

    // convert trajectories into robot frame
    cospsi = cos(psi_world);
    sinpsi = sin(psi_world);

    for(int i = 0; i < window_size; i++) 
    {
        x_veh[i] = traj_x[i] * cospsi + traj_y[i] * sinpsi - x_world;
        y_veh[i] = -traj_x[i] * sinpsi + traj_y[i] * cospsi - y_world;
    }

    // for(int i = 0; i < window_size; i++) 
    // {
    //     const double dx = traj_x[i] - cur_x;
    //     const double dy = traj_y[i] - cur_y;
    //     x_veh[i] = dx * cospsi + dy * sinpsi;
    //     y_veh[i] = dy * cospsi - dx * sinpsi;
    // }

    // 1rd order polynomial
    coeffs = polyfit(x_veh, y_veh, 1);
    // std::cout << "i : " << i << " / coeffs : " << coeffs << std::endl;

    const double cte  = polyeval(coeffs, 0.0);
    const double epsi = atan(coeffs[1]);

    // update state
    state << 0.0, 0.0, 0.0, cur_v, cur_w, cte, epsi;
  
    // for (auto x : traj_x)
    //     std::cout << "traj_x : " << x << std::endl;
    // for (auto y : traj_y)
    //     std::cout << "traj_y : " << y << std::endl;
  }

  auto gt_traj = GetTrajPointsCirc(0, iters, iters);
  auto gt_x = std::get<0>(gt_traj);
  auto gt_y = std::get<1>(gt_traj);

  if (true){
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

    plt::show();
  }

  return 0;
}