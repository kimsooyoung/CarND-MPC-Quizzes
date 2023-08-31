#include <matplotlibcpp.h>
#include <Eigen/QR>
#include <vector>
#include "circular_traj.h"
#include "helpers.h"
// #include "MPC.h"

namespace plt = matplotlibcpp;

using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

// For every steps
// 1. get traj
// 2. convert traj into robot frame
// 3. get coeffs
// 4. solve mpc
// 5. update states
// 6. plot results

int main() {
  // MPC mpc;

  const int window_size = 6;
  const int iters = 50;

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

  VectorXd state(6);
  state << x, y, psi, v, w, cte, epsi;

  for (size_t i = 0; i < iters; ++i) {
    cout << "Iteration " << i << endl;

    // auto vars = mpc.Solve(state, coeffs);
    auto vars = vector<double>(9, 0.0);

    cout << "x = " << vars[0] << endl;
    cout << "y = " << vars[1] << endl;
    cout << "psi = " << vars[2] << endl;
    cout << "v = " << vars[3] << endl;
    cout << "w = " << vars[4] << endl;
    cout << "cte = " << vars[5] << endl;
    cout << "epsi = " << vars[6] << endl;
    cout << "delta = " << vars[7] << endl;
    cout << "a = " << vars[8] << endl;
    cout << endl;

    auto cur_x = vars[0];
    auto cur_y = vars[1];
    auto cur_psi = vars[2];

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

  // // NOTE: free feel to play around with these
  // double x = -1;
  // double y = 10;
  // double psi = 0;
  // double v = 10;
  // // The cross track error is calculated by evaluating at polynomial at x, f(x)
  // // and subtracting y.
  // double cte = polyeval(coeffs, x) - y;
  // // Due to the sign starting at 0, the orientation error is -f'(x).
  // // derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
  // double epsi = psi - atan(coeffs[1]);

  // VectorXd state(6);
  // state << x, y, psi, v, cte, epsi;

  // vector<double> x_vals = {state[0]};
  // vector<double> y_vals = {state[1]};
  // vector<double> psi_vals = {state[2]};
  // vector<double> v_vals = {state[3]};
  // vector<double> cte_vals = {state[4]};
  // vector<double> epsi_vals = {state[5]};
  // vector<double> delta_vals = {};
  // vector<double> a_vals = {};

  // for (size_t i = 0; i < iters; ++i) {
  //   cout << "Iteration " << i << endl;

  //   auto vars = mpc.Solve(state, coeffs);

  //   x_vals.push_back(vars[0]);
  //   y_vals.push_back(vars[1]);
  //   psi_vals.push_back(vars[2]);
  //   v_vals.push_back(vars[3]);
  //   cte_vals.push_back(vars[4]);
  //   epsi_vals.push_back(vars[5]);

  //   delta_vals.push_back(vars[6]);
  //   a_vals.push_back(vars[7]);

  //   state << vars[0], vars[1], vars[2], vars[3], vars[4], vars[5];
  //   cout << "x = " << vars[0] << endl;
  //   cout << "y = " << vars[1] << endl;
  //   cout << "psi = " << vars[2] << endl;
  //   cout << "v = " << vars[3] << endl;
  //   cout << "cte = " << vars[4] << endl;
  //   cout << "epsi = " << vars[5] << endl;
  //   cout << "delta = " << vars[6] << endl;
  //   cout << "a = " << vars[7] << endl;
  //   cout << endl;
  // }

  // // TODO: matplotlibcpp  https://statphys.pknu.ac.kr/dokuwiki/doku.php?id=c:c_%EC%97%90_matplotlib_%EB%9D%BC%EC%9D%B4%EB%B8%8C%EB%9F%AC%EB%A6%AC_%EC%B6%94%EA%B0%80%ED%95%B4%EC%84%9C_%EA%B7%B8%EB%9E%98%ED%94%84_%EA%B7%B8%EB%A6%AC%EA%B8%B0
  // // Plot values
  // // NOTE: feel free to play around with this.
  // // It's useful for debugging!
  // plt::figure(1);
  // plt::subplot(3, 1, 1);
  // plt::title("X Values");
  // plt::plot(x_vals);
  // plt::subplot(3, 1, 2);
  // plt::title("Y Values");
  // plt::plot(y_vals);
  // plt::subplot(3, 1, 3);
  // plt::title("PSI Values");
  // plt::plot(psi_vals);

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

  // plt::show();

  return 0;
}