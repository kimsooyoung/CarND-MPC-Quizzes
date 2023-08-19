// In this quiz you'll fit a polynomial to waypoints.

#include <iostream>
#include <matplotlibcpp.h>
#include <Eigen/Core>
#include <Eigen/QR>

using Eigen::VectorXd;

namespace plt = matplotlibcpp;

// Evaluate a polynomial.
double polyeval(const VectorXd &coeffs, double x);
// Fit a polynomial.
VectorXd polyfit(const VectorXd &xvals, const VectorXd &yvals, int order);

int main() {

  std::vector<double> xvals_plot;
  std::vector<double> yvals_plot;

  // y = 2x + 1
  // VectorXd coeffs(2);
  // coeffs << 1, 2;

  // y = x(x-15)^2 / 20
  // y = x^3/20 - 30x^2/20 + 225x/20 + 0
  VectorXd coeffs(4);
  coeffs << 0.0, 225.0/20, -30.0/20, 1.0/20;

  for (double x = 0; x <= 20; ++x) {
    xvals_plot.push_back(x);
    yvals_plot.push_back(polyeval(coeffs, x));
  }

  plt::plot(xvals_plot, yvals_plot); //plot the x,y
  plt::grid(true); //show grid
  plt::show(); // show figure
}

double polyeval(const VectorXd &coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); ++i) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Adapted from:
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
VectorXd polyfit(const VectorXd &xvals, const VectorXd &yvals, int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);

  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); ++i) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); ++j) {
    for (int i = 0; i < order; ++i) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  std::cout << "A" << std::endl;
  std::cout << A << std::endl;

  std::cout << "yvals" << std::endl;
  std::cout << yvals << std::endl;

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);

  return result;
}