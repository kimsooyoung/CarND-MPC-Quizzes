#include "MPC.h"
#include <math.h>
#include <cppad/cppad.hpp>

#define HAVE_CSTDDEF
#include <cppad/ipopt/solve.hpp>
#undef HAVE_CSTDDEF

#include <vector>
#include <Eigen/Core>

using CppAD::AD;
using Eigen::VectorXd;

/**
 * TODO: Set N and dt
 */
size_t N = 5;
AD<double> dt = 0.2;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// NOTE: feel free to play around with this or do something completely different
double ref_v = 40;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  VectorXd coeffs;
  // Coefficients of the fitted polynomial.
  FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  // `fg` is a vector containing the cost and constraints.
  // `vars` is a vector containing the variable values (state & actuators).
  void operator()(ADvector& fg, const ADvector& vars) {
    assert( fg.size() == 1 + N * 6);
    assert( vars.size() == N * 6 + (N - 1) * 2 + 1);

    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    /*
     * TODO: Define the cost related the reference state and
     *   anything you think may be beneficial.
     */

    // for (size_t t = 0; t < N; ++t) {
    //   fg[0] += CppAD::pow(vars[cte_start + t], 2);
    //   fg[0] += CppAD::pow(vars[epsi_start + t], 2);
    //   fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    // }

    // for (size_t i = 0; i < N - 1; i++){
    //   fg[0] += CppAD::pow(vars[delta_start + i], 2);
    //   fg[0] += CppAD::pow(vars[a_start + i], 2);
    // }

    // for (size_t i = 0; i < N - 2; i++){
    //   fg[0] += CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
    //   fg[0] += CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    // }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    auto current_time = vars[N * 6 + (N - 1) * 2];

    // The rest of the constraints
    for (size_t t = 1; t < N; ++t) {
            
      /**
       * TODO: Grab the rest of the states at t+1 and t.
       *   We have given you parts of these states below.
       */
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // CppAD can compute derivatives and pass these to the solver.

      /**
       * TODO: Setup the rest of the model constraints
       */
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      
      // 1st order tajectory
      // f(x) = a0 + a1 * x
      //
      // 3rd order trajectory
      // f(x) = a0 + a1 * x + a2 * x^2 + a3 * x^3
      // f'(x) = a1 + 2 * a2 * x + 3 * a3 * x^2
      AD<double> psi_des = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0);
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psi_des) + v0 * delta0 / Lf * dt); 
      
      AD<double> f_x = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
      fg[1 + cte_start + t] = cte1 - (f_x - y0 + v0 * CppAD::sin(epsi0) * dt);
      
      // traj following mpc start
      std::cout << "current_time " << current_time << std::endl;

      // x_ref
      fg[0] += CppAD::pow(x1 - current_time, 2);
      // y_ref
      AD<double> f_y = coeffs[0] + coeffs[1] * x1 + coeffs[2] * x1 * x1 + coeffs[3] * x1 * x1 * x1;
      fg[0] += CppAD::pow(y1 - f_y, 2);
      // theta_ref
      // fg[0] += CppAD::pow(psi1, 2);

      // fg[1 + 6*N] = x1 - current_time;
      // fg[1 + 6*N + (N - 1)] = y1 - f_y;
      // fg[1 + 6*N + 2 * (N - 1)] = psi1 - 0;

      current_time += dt;
    }

    // std::cout << "fg.size() " << fg.size() << std::endl;
    // for(auto i : fg){
    //   std::cout << i << std::endl;
    // }
    // std::cout << "done" << std::endl;
  }
};

//
// MPC class definition
//

MPC::MPC() {}
MPC::~MPC() {}

std::vector<double> MPC::Solve(const VectorXd &x0, const VectorXd &coeffs) {
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x = x0[0];
  double y = x0[1];
  double psi = x0[2];
  double v = x0[3];
  double cte = x0[4];
  double epsi = x0[5];
  double t = x0[6];

  // number of independent variables
  // N timesteps == N - 1 actuations
  // time t added for trajectory
  size_t n_vars = N * 6 + (N - 1) * 2 + 1;
  // Number of constraints
  // trajectory added
  size_t n_constraints = N * 6;
  // size_t n_constraints = N * 6 + (N - 1) * 3;

  // Initial value of the independent variables.
  // Should be 0 except for the initial values.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; ++i) {
    vars[i] = 0.0;
  }

  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;
  vars[n_vars - 1] = t; 

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < delta_start; ++i) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (size_t i = delta_start; i < a_start; ++i) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (size_t i = a_start; i < n_vars; ++i) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  vars_lowerbound[n_vars - 1] = t; 
  vars_upperbound[n_vars - 1] = t; 

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; ++i) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  // x_0 < x0 < x_0
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // options
  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;


  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  //
  // Check some of the solution values
  //
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  return {solution.x[x_start + 1],   solution.x[y_start + 1],
          solution.x[psi_start + 1], solution.x[v_start + 1],
          solution.x[cte_start + 1], solution.x[epsi_start + 1],
          solution.x[delta_start],   solution.x[a_start]};
}