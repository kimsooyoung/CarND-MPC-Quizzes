#include <vector>
#include <matplotlibcpp.h>
#include <Eigen/QR>
#include "helpers.h"
#include "MPC.h"

namespace plt = matplotlibcpp;

using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

int main() {

  MPC mpc;
  // int iters = 50;
  int iters = 20;

  // multi point passing trajectory 
  VectorXd xvals(5);
  VectorXd yvals(5);

  // (0, 0) (5, 25) (10, 12.51) (15, 0) (20, 25)
  // y = a0 + a1x + a2x^2 + a3x^3
  xvals << 0.0, 5.0, 10.0, 15.0, 20.0;
  yvals << 0.0, 25.0, 12.51, 0.0, 25.0;
  auto coeffs = polyfit(xvals, yvals, 3);

  // initial state setup
  // NOTE: free feel to play around with these
  double x = -1;
  double y = 10;
  double psi = 0;
  double v = 10;
  /**
   * TODO: calculate the cross track error
   */
  double cte = polyeval(coeffs, x) - y;
  /**
   * TODO: calculate the orientation error
   */
  double epsi = psi - atan(coeffs[1]);

  VectorXd state(6);
  state << x, y, psi, v, cte, epsi;

  vector<double> x_vals = {state[0]};
  vector<double> y_vals = {state[1]};
  vector<double> psi_vals = {state[2]};
  vector<double> v_vals = {state[3]};
  vector<double> cte_vals = {state[4]};
  vector<double> epsi_vals = {state[5]};
  vector<double> delta_vals = {};
  vector<double> a_vals = {};

  for (int i = 0; i < iters; ++i) {
    // cout << "Iteration " << i << endl;

    auto vars = mpc.Solve(state, coeffs);

    x_vals.push_back(vars[0]);
    y_vals.push_back(vars[1]);
    psi_vals.push_back(vars[2]);
    v_vals.push_back(vars[3]);
    cte_vals.push_back(vars[4]);
    epsi_vals.push_back(vars[5]);

    delta_vals.push_back(vars[6]);
    a_vals.push_back(vars[7]);

    state << vars[0], vars[1], vars[2], vars[3], vars[4], vars[5];
    // cout << "x = " << vars[0] << endl;
    // cout << "y = " << vars[1] << endl;
    // cout << "psi = " << vars[2] << endl;
    // cout << "v = " << vars[3] << endl;
    // cout << "cte = " << vars[4] << endl;
    // cout << "epsi = " << vars[5] << endl;
    // cout << "delta = " << vars[6] << endl;
    // cout << "a = " << vars[7] << endl;
    // cout << endl;
  }

  // TODO: matplotlibcpp  https://statphys.pknu.ac.kr/dokuwiki/doku.php?id=c:c_%EC%97%90_matplotlib_%EB%9D%BC%EC%9D%B4%EB%B8%8C%EB%9F%AC%EB%A6%AC_%EC%B6%94%EA%B0%80%ED%95%B4%EC%84%9C_%EA%B7%B8%EB%9E%98%ED%94%84_%EA%B7%B8%EB%A6%AC%EA%B8%B0
  // Plot values
  // NOTE: feel free to play around with this.
  // It's useful for debugging!
  plt::figure(1);
  plt::subplot(3, 1, 1);
  // plt::title("XY Values");
  plt::plot(x_vals);
  plt::subplot(3, 1, 2);
  plt::plot(y_vals);
  // plt::title("PSI Values");
  plt::subplot(3, 1, 3);
  plt::plot(psi_vals);

  plt::figure(2);
  plt::subplot(3, 1, 1);
  plt::title("CTE");
  plt::plot(cte_vals);
  plt::subplot(3, 1, 2);
  plt::title("Delta (Radians)");
  plt::plot(delta_vals);
  plt::subplot(3, 1, 3);
  plt::title("Accel m/s^2");
  plt::plot(a_vals);

  plt::show();
}