// In this quiz you'll implement the global kinematic model.
#include <math.h>
#include <iostream>
#include "Eigen-3.3/Eigen/Core"

using Eigen::VectorXd;

//
// Helper functions
//
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const double Lf = 2;

// Return the next state.
VectorXd globalKinematic(const VectorXd &state, 
                         const VectorXd &actuators, double dt);

int main() {
  // [x, y, psi, v]
  VectorXd state(4);
  // [delta, v]
  VectorXd actuators(2);

  state << 0, 0, deg2rad(45), 1;
  actuators << deg2rad(5), 1;

  // should be [0.212132, 0.212132, 0.798488, 1.3]
  auto next_state = globalKinematic(state, actuators, 0.3);

  std::cout << next_state << std::endl;
}

VectorXd globalKinematic(const VectorXd &state, 
                         const VectorXd &actuators, double dt) {
  // Create a new vector for the next state.
  VectorXd next_state(state.size());

  /**
   * TODO: Implement the global kinematic model,
   *   to return the next state from the inputs.
   */
  double x_prev = state(0);
  double y_prev = state(1);
  double psi_prev = state(2);
  double v_prev = state(3);

  double delta_prev = actuators(0);
  double a_prev = actuators(1);

  double x_now = x_prev + v_prev * cos(psi_prev) * dt;
  double y_now = y_prev + v_prev * sin(psi_prev) * dt;
  double psi_now = psi_prev + v_prev / Lf * delta_prev * dt;
  double v_now = v_prev + a_prev * dt;

  // NOTE: state is [x, y, psi, v] and actuators is [delta, a]
  next_state << x_now, y_now, psi_now, v_now;

  return next_state;
}