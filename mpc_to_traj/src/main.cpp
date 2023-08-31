#include <matplotlibcpp.h>
#include <Eigen/QR>
#include <vector>
#include "circular_traj.h"
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

  const int window_size = 6;
  const int iters = 100;
  const double dt = 0.1;

  // initial trajectory
  auto traj_points = GetTrajPoints(0, window_size, iters);
  auto traj_x = std::get<0>(traj_points);
  auto traj_y = std::get<1>(traj_points);

  VectorXd ptsx(window_size);
  VectorXd ptsy(window_size);

  for(long unsigned int i = 0; i < window_size; ++i)
  {
    ptsx[i] = traj_x[i];
    ptsy[i] = traj_y[i];
  }

  // 3rd order polynomial
  auto coeffs = polyfit(ptsx, ptsy, 3);

  // initial state
  double x = 0.0;
  double y = 0.0;
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

  for (size_t i = 0; i < iters; ++i) {
    cout << "Iteration " << i << endl;

    auto vars = mpc.Solve(state, coeffs);

    // cout << "x = " << vars[0] << endl;
    // cout << "y = " << vars[1] << endl;
    // cout << "psi = " << vars[2] << endl;
    // cout << "v = " << vars[3] << endl;
    // cout << "w = " << vars[4] << endl;
    // cout << "cte = " << vars[5] << endl;
    // cout << "epsi = " << vars[6] << endl;
    // cout << "delta = " << vars[7] << endl;
    // cout << "a = " << vars[8] << endl;
    // cout << endl;

    auto cur_x = vars[0];
    auto cur_y = vars[1];
    auto cur_psi = vars[2];
    auto cur_v = vars[3];
    auto cur_w = vars[4];

    x_vals.push_back(cur_x + cur_v * cos(cur_psi) * dt);
    y_vals.push_back(cur_y + cur_v * sin(cur_psi) * dt);
    psi_vals.push_back(cur_psi + cur_w * dt);

    // x_vals.push_back(cur_x);
    // y_vals.push_back(cur_y);
    // psi_vals.push_back(cur_psi + cur_w * dt);

    // update trajectory
    auto traj_points = GetTrajPoints(i, window_size, iters);
    auto traj_x = std::get<0>(traj_points);
    auto traj_y = std::get<1>(traj_points);

    // Convert to the vehicle coordinate system
    VectorXd x_veh(window_size);
    VectorXd y_veh(window_size);

    const double cospsi = cos(cur_psi);
    const double sinpsi = sin(cur_psi);

    for(int i = 0; i < window_size; i++) 
    {
        const double dx = traj_x[i] - cur_x;
        const double dy = traj_x[i] - cur_y;
        x_veh[i] = dx * cospsi + dy * sinpsi;
        y_veh[i] = dy * cospsi - dx * sinpsi;
    }
    
    // Fit waypoints
    coeffs = polyfit(x_veh, y_veh, 3); 

    const double cte  = polyeval(coeffs, 0.0);
    const double epsi = atan(coeffs[1]);

    // update state
    state << vars[0], vars[1], vars[2], vars[3], vars[4], cte, epsi;
  }

  const double pi = M_PI;
  vector<double> gt_x(50 + 1);
  vector<double> gt_y(50 + 1);

  for(int i = 0; i < 50 + 1; ++i)
  {
      gt_x[i] = cos(2 * pi * i / 50);
      gt_y[i] = sin(2 * pi * i / 50) + 1;
  }

  // plt::plot(gt_x, gt_y, "r--"); //plot the x,y
  // plt::plot(x_vals, y_vals); //plot the x,y
  // plt::grid(true); //show grid
  // plt::show(); // show figure

  // TODO: matplotlibcpp  https://statphys.pknu.ac.kr/dokuwiki/doku.php?id=c:c_%EC%97%90_matplotlib_%EB%9D%BC%EC%9D%B4%EB%B8%8C%EB%9F%AC%EB%A6%AC_%EC%B6%94%EA%B0%80%ED%95%B4%EC%84%9C_%EA%B7%B8%EB%9E%98%ED%94%84_%EA%B7%B8%EB%A6%AC%EA%B8%B0
  // Plot values
  // NOTE: feel free to play around with this.
  // It's useful for debugging!
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

  // plt::figure(2);
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

  return 0;
}