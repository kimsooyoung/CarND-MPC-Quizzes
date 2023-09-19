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

class FG_eval {
private:
  int mpc_step_;
  int x_start_, y_start_, psi_start_;
  int v_start_, w_start_; 
  int cte_start_, epsi_start_;
  int a_start_, alpha_start_;

  double _w_cte, _w_epsi, _w_vel;
  double _w_a, _w_alpha;
  double _w_delta_a, _w_delta_alpha;

  double dt;
  double ref_v_;

public:
  VectorXd coeffs_;

  FG_eval(VectorXd coeffs) { 
    coeffs_ = coeffs; 
    
    mpc_step_ = 40;

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
    mpc_step_ = params.find("STEPS") != params.end() ? params.at("STEPS") : mpc_step_;
    ref_v_   = params.find("REF_V") != params.end() ? params.at("REF_V") : ref_v_;
    dt = params.find("DT") != params.end() ? params.at("DT") : dt;
    
    _w_cte   = params.find("W_CTE") != params.end()   ? params.at("W_CTE") : _w_cte;
    _w_epsi  = params.find("W_EPSI") != params.end()  ? params.at("W_EPSI") : _w_epsi;
    _w_vel   = params.find("W_V") != params.end()     ? params.at("W_V") : _w_vel;

    _w_a = params.find("W_A") != params.end()     ? params.at("W_A") : _w_a;
    _w_alpha = params.find("W_ALPHA") != params.end() ? params.at("W_ALPHA") : _w_alpha;
    _w_delta_a = params.find("W_DELTA_A") != params.end() ? params.at("W_DELTA_A") : _w_delta_a;
    _w_delta_alpha = params.find("W_DELTA_ALPHA") != params.end() ? params.at("W_DELTA_ALPHA") : _w_delta_alpha;

    std::cout << "[FG_eval] mpc_step : " << mpc_step_ << std::endl;
    std::cout << "[FG_eval] ref_v : " << ref_v_ << std::endl;
    std::cout << "[FG_eval] dt : " << dt << std::endl;

    std::cout << "[FG_eval] _w_cte : " << _w_cte << std::endl;
    std::cout << "[FG_eval] _w_epsi : " << _w_epsi << std::endl;
    std::cout << "[FG_eval] _w_vel : " << _w_vel << std::endl;

    std::cout << "[FG_eval] _w_a : " << _w_a << std::endl;
    std::cout << "[FG_eval] _w_alpha : " << _w_alpha << std::endl;
    std::cout << "[FG_eval] _w_delta_a : " << _w_delta_a << std::endl;
    std::cout << "[FG_eval] _w_delta_alpha : " << _w_delta_alpha << std::endl;

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

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    fg[0] = 0;
    
    for (auto i = 0; i < mpc_step_; ++i) {
      fg[0] += _w_cte * CppAD::pow(vars[cte_start_ + i], 2);
      fg[0] += _w_epsi * CppAD::pow(vars[epsi_start_ + i], 2);
      fg[0] += _w_vel * CppAD::pow(vars[v_start_ + i] - ref_v_, 2);
    }

    for (auto i = 0; i < mpc_step_ - 1; i++){
      fg[0] += _w_a * CppAD::pow(vars[a_start_ + i], 2);
      fg[0] += _w_alpha * CppAD::pow(vars[alpha_start_ + i], 2);
    }

    for (auto i = 0; i < mpc_step_ - 2; i++){
      fg[0] += _w_delta_a * CppAD::pow(vars[a_start_ + i + 1] - vars[a_start_ + i], 2);
      fg[0] += _w_delta_alpha * CppAD::pow(vars[alpha_start_ + i + 1] - vars[alpha_start_ + i], 2);
    }

    fg[1 + x_start_] = vars[x_start_];
    fg[1 + y_start_] = vars[y_start_];
    fg[1 + psi_start_] = vars[psi_start_];
    fg[1 + v_start_] = vars[v_start_];
    fg[1 + w_start_] = vars[w_start_];
    fg[1 + cte_start_] = vars[cte_start_];
    fg[1 + epsi_start_] = vars[epsi_start_];

    for (auto t = 1; t < mpc_step_; ++t) {

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
      for (int i = 0; i < coeffs_.size(); i++)
        f0 += coeffs_[i] * CppAD::pow(x0, i);

      AD<double> psides0 = 0.0;
      for (int i = 1; i < coeffs_.size(); i++) 
          psides0 += i * coeffs_[i] * CppAD::pow(x0, i-1); // f'(x0)

      if (abs(psides0) < 1e-3)
        psides0 = 0.0;
      else
        psides0 = CppAD::atan(psides0);

      fg[1 + x_start_ + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start_ + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start_ + t] = psi1 - (psi0 + w0 * dt);
      fg[1 + v_start_ + t] = v1 - (v0 + a0 * dt);
      fg[1 + w_start_ + t] = w1 - (w0 + alpha0 * dt);
      fg[1 + cte_start_ + t] = cte1 - ((f0 - y0));
      fg[1 + epsi_start_ + t] = epsi1 - ((psi0 - psides0));
    }
  }
};

MPC::MPC() {}
MPC::~MPC() {}

void MPC::LoadParams(const std::map<std::string, double> &params){

  params_ = params;

  mpc_step_ = params.find("STEPS") != params.end() ? params.at("STEPS") : mpc_step_;
  min_acc_ = params.find("MIN_ACC") != params.end() ? params.at("MIN_ACC") : mpc_step_;
  max_acc_ = params.find("MAX_ACC") != params.end() ? params.at("MAX_ACC") : mpc_step_;
  min_ang_acc_ = params.find("MIN_ANG_ACC") != params.end() ? params.at("MIN_ANG_ACC") : mpc_step_;
  max_ang_acc_ = params.find("MAX_ANG_ACC") != params.end() ? params.at("MAX_ANG_ACC") : mpc_step_;
  
  x_start_     = 0;
  y_start_     = x_start_ + mpc_step_;
  psi_start_   = y_start_ + mpc_step_;
  v_start_     = psi_start_ + mpc_step_;
  w_start_     = v_start_ + mpc_step_;
  cte_start_   = w_start_ + mpc_step_;
  epsi_start_  = cte_start_ + mpc_step_;
  a_start_     = epsi_start_ + mpc_step_;
  alpha_start_ = a_start_ + mpc_step_ - 1;

  std::cout << "[MPC] mpc_step : " << mpc_step_ << std::endl;
}


std::vector<double> MPC::Solve(const VectorXd &x0, const VectorXd &coeffs) {

  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  double x = x0[0];
  double y = x0[1];
  double psi = x0[2];
  double v_x = x0[3];
  double v_y = x0[4];
  double w = x0[5];
  double cte_x = x0[6];
  double cte_y = x0[7];
  
  int n_vars = mpc_step_ * 8 + (mpc_step_ - 1) * 2;
  Dvector vars(n_vars);

  for (auto i = 0; i < n_vars; ++i) {
    vars[i] = 0.0;
  }

  vars[x_start_] = x;
  vars[y_start_] = y;
  vars[psi_start_] = psi;
  vars[v_start_] = v;
  vars[w_start_] = w;
  vars[cte_start_] = cte;
  vars[epsi_start_] = epsi;

  Dvector vars_lowerbound(n_vars), vars_upperbound(n_vars);

  for (auto i = 0; i < a_start_; ++i) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (auto i = a_start_; i < alpha_start_; ++i) {
    vars_lowerbound[i] = min_acc_;
    vars_upperbound[i] = max_acc_;
  }

  for (auto i = alpha_start_; i < n_vars; ++i) {
    vars_lowerbound[i] = min_ang_acc_;
    vars_upperbound[i] = max_ang_acc_;
  }

  auto n_constraints = mpc_step_ * 7;

  Dvector constraints_lowerbound(n_constraints), constraints_upperbound(n_constraints);

  for (auto i = 0; i < n_constraints; ++i) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

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

  FG_eval fg_eval(coeffs);
  fg_eval.LoadParams(params_);

  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";

  CppAD::ipopt::solve_result<Dvector> solution;
  
  CppAD::ipopt::solve<Dvector, FG_eval>(
    options, vars, vars_lowerbound, vars_upperbound, 
    constraints_lowerbound, constraints_upperbound, 
    fg_eval, solution
  );

  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  
  // auto cost = solution.obj_value;
  auto answer = solution.x;

  return {answer[x_start_ + 1],   answer[y_start_ + 1],
          answer[psi_start_ + 1], answer[v_start_ + 1],
          answer[w_start_ + 1],   answer[cte_start_ + 1],
          answer[epsi_start_ + 1],answer[a_start_],
          answer[alpha_start_]};
}