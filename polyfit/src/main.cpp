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
  VectorXd xvals(5);
  VectorXd yvals(5);

  // (0, 0) (5, 25) (10, 12.51) (15, 0) (20, 25)
  // y = a0 + a1x + a2x^2 + a3x^3
  // x waypoint coordinates
  xvals << 0.0, 5.0, 10.0, 15.0, 20.0;
  yvals << 0.0, 25.0, 12.51, 0.0, 25.0;

  /**
   * TODO: use `polyfit` to fit a third order polynomial to the (x, y)
   *   coordinates.
   * Hint: call Eigen::VectorXd polyfit() and pass xvals, yvals, and the 
   *   polynomial degree/order
   */
  // YOUR CODE HERE
  // VectorXd polyfit(const VectorXd &xvals, const VectorXd &yvals, int order)
  VectorXd coeffs = polyfit(xvals, yvals, 3);
  std::cout << "Coefficients: " << std::endl;
  std::cout << coeffs << std::endl << std::endl << std::endl;

  std::vector<double> xvals_plot;
  std::vector<double> yvals_plot;

  for (double x = 0; x <= 20; ++x) {
    /**
     * TODO: use `polyeval` to evaluate the x values.
     */
    xvals_plot.push_back(x);
    yvals_plot.push_back(polyeval(coeffs, x));
    // std::cout << polyeval(coeffs, x) << std::endl;
  }

  plt::plot(xvals_plot, yvals_plot); //plot the x,y
  plt::grid(true); //show grid
  plt::show(); // show figure
  
  // Expected output
  // -0.905562
  // -0.226606
  // 0.447594
  // 1.11706
  // 1.7818
  // 2.44185
  // 3.09723
  // 3.74794
  // 4.39402
  // 5.03548
  // 5.67235
  // 6.30463
  // 6.93236
  // 7.55555
  // 8.17423
  // 8.7884
  // 9.3981
  // 10.0033
  // 10.6041
  // 11.2005
  // 11.7925
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
  std::cout << A << std::endl << std::endl;

  std::cout << "yvals" << std::endl;
  std::cout << yvals << std::endl << std::endl;

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);

  return result;
}