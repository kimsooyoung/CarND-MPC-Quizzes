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
size_t N = 25;
double dt = 0.05;

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
// size_t x_start = 0;
// size_t y_start = x_start + N;
// size_t psi_start = y_start + N;
// size_t v_start = psi_start + N;
// size_t cte_start = v_start + N;
// size_t epsi_start = cte_start + N;
// size_t delta_start = epsi_start + N;
// size_t a_start = delta_start + N - 1;

class FG_eval {
private:
  int mpc_step_;
  int x_start_, y_start_, psi_start_;
  int v_start_, w_start_; 
  int cte_start_, epsi_start_;
  int a_start_, alpha_start_;

  double dt;
  double ref_v_;

public:
  VectorXd coeffs;
  // Coefficients of the fitted polynomial.
  FG_eval(VectorXd coeffs) { 
    this->coeffs = coeffs; 

    dt = 0.1;
    mpc_step_ = 40;
    ref_v_ = 0.5;

    x_start_     = 0;
    y_start_     = x_start_ + mpc_step_;
    psi_start_   = y_start_ + mpc_step_;
    v_start_     = psi_start_ + mpc_step_;
    w_start_     = v_start_ + mpc_step_;
    cte_start_   = w_start_ + mpc_step_;
    epsi_start_  = cte_start_ + mpc_step_;
    a_start_     = epsi_start_ + mpc_step_;
    alpha_start_ = a_start_ + mpc_step_ - 1;
  }

  void LoadParams(const std::map<std::string, double> &params){
    dt = params.find("DT") != params.end() ? params.at("DT") : dt;
    mpc_step_ = params.find("STEPS") != params.end() ? params.at("STEPS") : mpc_step_;
    ref_v_   = params.find("REF_V") != params.end() ? params.at("REF_V") : ref_v_;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  // `fg` is a vector containing the cost and constraints.
  // `vars` is a vector containing the variable values (state & actuators).
  void operator()(ADvector& fg, const ADvector& vars) {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    /*
     * TODO: Define the cost related the reference state and
     *   anything you think may be beneficial.
     */

    for (auto t = 0; t < mpc_step_; ++t) {
      fg[0] += CppAD::pow(vars[cte_start_ + t], 2);
      fg[0] += CppAD::pow(vars[epsi_start_ + t], 2);
      fg[0] += CppAD::pow(vars[v_start_ + t] - ref_v_, 2);
    }

    for (auto i = 0; i < mpc_step_ - 1; i++){
      fg[0] += CppAD::pow(vars[a_start_ + i], 2);
      fg[0] += CppAD::pow(vars[alpha_start_ + i], 2);
    }

    for (auto i = 0; i < mpc_step_ - 2; i++){
      fg[0] += CppAD::pow(vars[a_start_ + i + 1] - vars[a_start_ + i], 2);
      fg[0] += CppAD::pow(vars[alpha_start_ + i + 1] - vars[alpha_start_ + i], 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start_] = vars[x_start_];
    fg[1 + y_start_] = vars[y_start_];
    fg[1 + psi_start_] = vars[psi_start_];
    fg[1 + v_start_] = vars[v_start_];
    fg[1 + w_start_] = vars[w_start_];
    fg[1 + cte_start_] = vars[cte_start_];
    fg[1 + epsi_start_] = vars[epsi_start_];

    // The rest of the constraints
    for (auto t = 1; t < mpc_step_; ++t) {
      /**
       * TODO: Grab the rest of the states at t+1 and t.
       *   We have given you parts of these states below.
       */
      AD<double> x1 = vars[x_start_ + t];
      AD<double> y1 = vars[y_start_ + t];
      AD<double> psi1 = vars[psi_start_ + t];
      AD<double> v1 = vars[v_start_ + t];
      AD<double> w1 = vars[w_start_ + t];
      AD<double> cte1 = vars[cte_start_ + t];
      AD<double> epsi1 = vars[epsi_start_ + t];

      AD<double> x0 = vars[x_start_ + t - 1];
      AD<double> y0 = vars[y_start_ + t - 1];
      AD<double> psi0 = vars[psi_start_ + t - 1];
      AD<double> v0 = vars[v_start_ + t - 1];
      AD<double> w0 = vars[w_start_ + t - 1];
      AD<double> cte0 = vars[cte_start_ + t - 1];
      AD<double> epsi0 = vars[epsi_start_ + t - 1];

      AD<double> a0 = vars[a_start_ + t - 1];
      AD<double> alpha0 = vars[alpha_start_ + t - 1];

      AD<double> f0 = 0.0;
      for (int i = 0; i < coeffs.size(); i++) 
          f0 += coeffs[i] * CppAD::pow(x0, i);

      AD<double> psides0 = 0.0;
      for (int i = 1; i < coeffs.size(); i++) 
          psides0 += i * coeffs[i] * CppAD::pow(x0, i-1); // f'(x0)
      psides0 = CppAD::atan(psides0);

      fg[1 + x_start_ + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start_ + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start_ + t] = psi1 - (psi0 + w0 * dt);
      fg[1 + v_start_ + t] = v1 - (v0 + a0 * dt);
      fg[1 + w_start_ + t] = w1 - (w0 + alpha0 * dt);
      // fg[1 + cte_start + t] =
      //     cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      // fg[1 + epsi_start + t] =
      //     epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);

      fg[1 + cte_start_ + t] = cte1 - (f0 - y0);
      fg[1 + epsi_start_ + t] = epsi1 - (psi0 - psides0);
    }
  }
};

//
// MPC class definition
//

MPC::MPC() {}
MPC::~MPC() {}

void MPC::LoadParams(const std::map<std::string, double> &params){

  params_ = params;
  mpc_step_ = params.find("STEPS") != params.end() ? params.at("STEPS") : mpc_step_;

  x_start_     = 0;
  y_start_     = x_start_ + mpc_step_;
  psi_start_   = y_start_ + mpc_step_;
  v_start_     = psi_start_ + mpc_step_;
  w_start_     = v_start_ + mpc_step_;
  cte_start_   = w_start_ + mpc_step_;
  epsi_start_  = cte_start_ + mpc_step_;
  a_start_     = epsi_start_ + mpc_step_;
  alpha_start_ = a_start_ + mpc_step_ - 1;
}


std::vector<double> MPC::Solve(const VectorXd &x0, const VectorXd &coeffs) {
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // set initial state 
  // set state ub, lb
  // set constraints ub, lb
  // solve
  // prepare return values

  double x = x0[0];
  double y = x0[1];
  double psi = x0[2];
  double v = x0[3];
  double w = x0[4];
  double cte = x0[5];
  double epsi = x0[6];

  // number of independent variables
  // N timesteps == N - 1 actuations
  size_t n_vars = mpc_step_ * 7 + (mpc_step_ - 1) * 2;

  // Initial value of the independent variables.
  // Should be 0 except for the initial values.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; ++i) {
    vars[i] = 0.0;
  }
  // Set the initial variable values
  vars[x_start_] = x;
  vars[y_start_] = y;
  vars[psi_start_] = psi;
  vars[v_start_] = v;
  vars[w_start_] = w;
  vars[cte_start_] = cte;
  vars[epsi_start_] = epsi;
  vars[a_start_] = 0.0;
  vars[alpha_start_] = 0.0;

  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (auto i = 0; i < a_start_; ++i) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  // min, max linear acc : 2.5 m/s^2
  for (auto i = a_start_; i < alpha_start_; ++i) {
    vars_lowerbound[i] = -2.5;
    vars_upperbound[i] = 2.5;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  // min, max angular acc : 1.2 rad/s^2
  for (size_t i = alpha_start_; i < n_vars; ++i) {
    vars_lowerbound[i] = -1.2;
    vars_upperbound[i] = 1.2;
  }

  // Number of constraints
  // 7 initial cond + dynamics conds.
  // eventually 7 * mpc_step_
  size_t n_constraints = 7 * (mpc_step_ - 1) +  7;

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
  constraints_lowerbound[x_start_] = x;
  constraints_lowerbound[y_start_] = y;
  constraints_lowerbound[psi_start_] = psi;
  constraints_lowerbound[v_start_] = v;
  constraints_lowerbound[w_start_] = w;
  constraints_lowerbound[cte_start_] = cte;
  constraints_lowerbound[epsi_start_] = epsi;

  constraints_upperbound[x_start_] = x;
  constraints_upperbound[y_start_] = y;
  constraints_upperbound[psi_start_] = psi;
  constraints_upperbound[v_start_] = v;
  constraints_upperbound[w_start_] = w;
  constraints_upperbound[cte_start_] = cte;
  constraints_upperbound[epsi_start_] = epsi;

  std::cout << "constraints_upperbound.size() : " << constraints_upperbound.size() << std::endl;

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs);
  fg_eval.LoadParams(params_);

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
  std::cout << solution.x.size() << std::endl;

  return {solution.x[x_start_ + 1],   solution.x[y_start_ + 1],
          solution.x[psi_start_ + 1], solution.x[v_start_ + 1],
          solution.x[w_start_ + 1], solution.x[cte_start_ + 1], 
          solution.x[epsi_start_ + 1], solution.x[a_start_],   
          solution.x[alpha_start_]};
}