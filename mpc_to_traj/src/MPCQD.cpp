#include "MPCQD.h"
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
  int vx_start_, vy_start_, w_start_; 
  int cte_x_start_, cte_y_start_;
  int ax_start_, ay_start_, alpha_start_;

  double _w_cte_x, _w_cte_y, _w_vel;
  double _w_ax, _w_ay, _w_alpha;
  double _w_delta_ax, _w_delta_ay, _w_delta_alpha;

  double dt;
  double ref_v_;

public:
  vector<double> traj_x_;
  vector<double> traj_y_;

  FG_eval(tuple<vector<double>, vector<double>> trajs) { 
    traj_x_ = std::get<0>(trajs);
    traj_y_ = std::get<1>(trajs);

    // for(auto x : traj_x_)
    //   std::cout << x << " ";
    // std::cout << std::endl;

    // for(auto y: traj_y_)
    //   std::cout << y << " ";
    // std::cout << std::endl;

    mpc_step_ = 40;

    x_start_     = 0;
    y_start_     = x_start_   + mpc_step_;
    psi_start_   = y_start_   + mpc_step_;
    vx_start_    = psi_start_ + mpc_step_;
    vy_start_    = vx_start_  + mpc_step_;
    w_start_     = vy_start_  + mpc_step_;
    cte_x_start_ = w_start_   + mpc_step_;
    cte_y_start_ = cte_x_start_ + mpc_step_;
    ax_start_    = cte_y_start_ + mpc_step_;
    ay_start_    = ax_start_    + mpc_step_ - 1;
    alpha_start_ = ay_start_    + mpc_step_ - 1;
  }

  void LoadParams(const std::map<std::string, double> &params){
    mpc_step_ = params.find("STEPS") != params.end() ? params.at("STEPS") : mpc_step_;
    ref_v_   = params.find("REF_V") != params.end() ? params.at("REF_V") : ref_v_;
    dt = params.find("DT") != params.end() ? params.at("DT") : dt;
    
    _w_cte_x   = params.find("W_CTE_X") != params.end() ? params.at("W_CTE_X") : _w_cte_x;
    _w_cte_y   = params.find("W_CTE_Y") != params.end() ? params.at("W_CTE_Y") : _w_cte_y;
    _w_vel   = params.find("W_V") != params.end() ? params.at("W_V") : _w_vel;

    _w_ax = params.find("W_AX") != params.end() ? params.at("W_AX") : _w_ax;
    _w_ay = params.find("W_AY") != params.end() ? params.at("W_AY") : _w_ay;
    _w_alpha = params.find("W_ALPHA") != params.end() ? params.at("W_ALPHA") : _w_alpha;
    
    _w_delta_ax = params.find("W_DELTA_AX") != params.end() ? params.at("W_DELTA_AX") : _w_delta_ax;
    _w_delta_ay = params.find("W_DELTA_AY") != params.end() ? params.at("W_DELTA_AY") : _w_delta_ay;
    _w_delta_alpha = params.find("W_DELTA_ALPHA") != params.end() ? params.at("W_DELTA_ALPHA") : _w_delta_alpha;

    // std::cout << "[FG_eval] mpc_step : " << mpc_step_ << std::endl;
    // std::cout << "[FG_eval] ref_v : " << ref_v_ << std::endl;
    // std::cout << "[FG_eval] dt : " << dt << std::endl;

    // std::cout << "[FG_eval] _w_cte_x : " << _w_cte_x << std::endl;
    // std::cout << "[FG_eval] _w_cte_y : " << _w_cte_y << std::endl;
    // std::cout << "[FG_eval] _w_vel : " << _w_vel << std::endl;

    // std::cout << "[FG_eval] _w_ax : " << _w_ax << std::endl;
    // std::cout << "[FG_eval] _w_ay : " << _w_ay << std::endl;
    // std::cout << "[FG_eval] _w_alpha : " << _w_alpha << std::endl;
    
    // std::cout << "[FG_eval] _w_delta_ax : " << _w_delta_ax << std::endl;
    // std::cout << "[FG_eval] _w_delta_ay : " << _w_delta_ay << std::endl;
    // std::cout << "[FG_eval] _w_delta_alpha : " << _w_delta_alpha << std::endl;

    x_start_     = 0;
    y_start_     = x_start_   + mpc_step_;
    psi_start_   = y_start_   + mpc_step_;
    vx_start_    = psi_start_ + mpc_step_;
    vy_start_    = vx_start_  + mpc_step_;
    w_start_     = vy_start_  + mpc_step_;
    cte_x_start_ = w_start_   + mpc_step_;
    cte_y_start_ = cte_x_start_ + mpc_step_;
    ax_start_    = cte_y_start_ + mpc_step_;
    ay_start_    = ax_start_    + mpc_step_ - 1;
    alpha_start_ = ay_start_    + mpc_step_ - 1;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    fg[0] = 0;
    
    for (auto i = 0; i < mpc_step_; ++i) {
      fg[0] += _w_cte_x * CppAD::pow(vars[cte_x_start_ + i], 2);
      fg[0] += _w_cte_y * CppAD::pow(vars[cte_y_start_ + i], 2);

      auto vx = vars[vx_start_ + i];
      auto vy = vars[vy_start_ + i];
      // CppAD::sqrt has nan error
      fg[0] += _w_vel * CppAD::pow( 
        (vx * vx + vy * vy) - CppAD::pow(ref_v_, 2), 2
      );
    }

    for (auto i = 0; i < mpc_step_ - 1; i++){
      fg[0] += _w_ax * CppAD::pow(vars[ax_start_ + i], 2);
      fg[0] += _w_ay * CppAD::pow(vars[ay_start_ + i], 2);
      fg[0] += _w_alpha * CppAD::pow(vars[alpha_start_ + i], 2);
    }

    for (auto i = 0; i < mpc_step_ - 2; i++){
      fg[0] += _w_delta_ax * CppAD::pow(vars[ax_start_ + i + 1] - vars[ax_start_ + i], 2);
      fg[0] += _w_delta_ay * CppAD::pow(vars[ay_start_ + i + 1] - vars[ay_start_ + i], 2);
      fg[0] += _w_delta_alpha * CppAD::pow(vars[alpha_start_ + i + 1] - vars[alpha_start_ + i], 2);
    }

    fg[1 + x_start_] = vars[x_start_];
    fg[1 + y_start_] = vars[y_start_];
    fg[1 + psi_start_] = vars[psi_start_];
    fg[1 + vx_start_] = vars[vx_start_];
    fg[1 + vy_start_] = vars[vy_start_];
    fg[1 + w_start_] = vars[w_start_];
    fg[1 + cte_x_start_] = vars[cte_x_start_];
    fg[1 + cte_y_start_] = vars[cte_y_start_];

    for (auto t = 1; t < mpc_step_; ++t) {

      AD<double> x1 = vars[x_start_ + t];
      AD<double> y1 = vars[y_start_ + t];
      AD<double> psi1 = vars[psi_start_ + t];
      AD<double> vx1 = vars[vx_start_ + t];
      AD<double> vy1 = vars[vy_start_ + t];
      AD<double> w1 = vars[w_start_ + t];
      AD<double> cte_x1 = vars[cte_x_start_ + t];
      AD<double> cte_y1 = vars[cte_y_start_ + t];

      AD<double> x0 = vars[x_start_ + t - 1];
      AD<double> y0 = vars[y_start_ + t - 1];
      AD<double> psi0 = vars[psi_start_ + t - 1];
      AD<double> vx0 = vars[vx_start_ + t - 1];
      AD<double> vy0 = vars[vy_start_ + t - 1];
      AD<double> w0 = vars[w_start_ + t - 1];
      AD<double> cte_x0 = vars[cte_x_start_ + t - 1];
      AD<double> cte_y0 = vars[cte_y_start_ + t - 1];

      AD<double> ax0 = vars[ax_start_ + t - 1];
      AD<double> ay0 = vars[ay_start_ + t - 1];
      AD<double> alpha0 = vars[alpha_start_ + t - 1];

      fg[1 + x_start_ + t] = x1 - (x0 + vx0 * dt);
      fg[1 + y_start_ + t] = y1 - (y0 + vy0 * dt);
      fg[1 + psi_start_ + t] = psi1 - (psi0 + w0 * dt);
      fg[1 + vx_start_ + t] = vx1 - (vx0 + ax0 * dt);
      fg[1 + vy_start_ + t] = vy1 - (vy0 + ay0 * dt);
      fg[1 + w_start_ + t] = w1 - (w0 + alpha0 * dt);
      fg[1 + cte_x_start_ + t] = cte_x1 - (traj_x_[t-1] - x0);
      fg[1 + cte_y_start_ + t] = cte_y1 - (traj_y_[t-1] - y0);

      // std::cout << "traj_x_[t-1] - cte_x0 : " << traj_x_[t-1] - cte_x0 << std::endl;
    }
  }
};

MPCQD::MPCQD() {}
MPCQD::~MPCQD() {}

void MPCQD::LoadParams(const std::map<std::string, double> &params){

  params_ = params;

  mpc_step_ = params.find("STEPS") != params.end() ? params.at("STEPS") : mpc_step_;
  
  min_acc_x_ = params.find("MIN_ACC_X") != params.end() ? params.at("MIN_ACC_X") : min_acc_x_;
  max_acc_x_ = params.find("MAX_ACC_X") != params.end() ? params.at("MAX_ACC_X") : max_acc_x_;
  
  min_acc_y_ = params.find("MIN_ACC_Y") != params.end() ? params.at("MIN_ACC_Y") : min_acc_x_;
  max_acc_y_ = params.find("MAX_ACC_Y") != params.end() ? params.at("MAX_ACC_Y") : max_acc_y_;
  
  min_ang_acc_ = params.find("MIN_ANG_ACC") != params.end() ? params.at("MIN_ANG_ACC") : min_ang_acc_;
  max_ang_acc_ = params.find("MAX_ANG_ACC") != params.end() ? params.at("MAX_ANG_ACC") : max_ang_acc_;
  
  x_start_     = 0;
  y_start_     = x_start_ + mpc_step_;
  psi_start_   = y_start_ + mpc_step_;
  vx_start_    = psi_start_ + mpc_step_;
  vy_start_    = vx_start_ + mpc_step_;
  w_start_     = vy_start_ + mpc_step_;
  cte_x_start_ = w_start_ + mpc_step_;
  cte_y_start_ = cte_x_start_ + mpc_step_;
  
  ax_start_    = cte_y_start_ + mpc_step_;
  ay_start_    = ax_start_ + mpc_step_ - 1;
  alpha_start_ = ay_start_ + mpc_step_ - 1;

  std::cout << "[MPC] mpc_step : " << mpc_step_ << std::endl;
  std::cout << "[MPC] min_acc_x : " << min_acc_x_ << std::endl;
  std::cout << "[MPC] max_acc_x : " << max_acc_x_ << std::endl;
  std::cout << "[MPC] min_acc_y : " << min_acc_y_ << std::endl;
  std::cout << "[MPC] max_acc_y : " << max_acc_y_ << std::endl;
  std::cout << "[MPC] min_ang_acc : " << min_ang_acc_ << std::endl;
  std::cout << "[MPC] max_ang_acc : " << max_ang_acc_ << std::endl;
}

std::vector<double> MPCQD::Solve(const VectorXd &x0, const tuple<vector<double>, vector<double>> &trajs) {

  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  // state vector dimension
  const int state_dim = 8;
  const int actuation_dim = 3;

  double x = x0[0];
  double y = x0[1];
  double psi = x0[2];
  double vx = x0[3];
  double vy = x0[4];
  double w = x0[5];
  double cte_x = x0[6];
  double cte_y = x0[7];
  
  // (8 state vars) * N + (3 actuator vars) * (N - 1) 
  int n_vars = mpc_step_ * state_dim + (mpc_step_ - 1) * actuation_dim;
  Dvector vars(n_vars);

  for (auto i = 0; i < n_vars; ++i) {
    vars[i] = 0.0;
  }

  vars[x_start_] = x;
  vars[y_start_] = y;
  vars[psi_start_] = psi;
  vars[vx_start_] = vx;
  vars[vy_start_] = vy;
  vars[w_start_] = w;
  vars[cte_x_start_] = cte_x;
  vars[cte_y_start_] = cte_y;

  // boundary conditions
  Dvector vars_lowerbound(n_vars), vars_upperbound(n_vars);

  for (auto i = 0; i < ax_start_; ++i) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (auto i = ax_start_; i < ay_start_; ++i) {
    vars_lowerbound[i] = min_acc_x_;
    vars_upperbound[i] = max_acc_x_;
  }

  for (auto i = ay_start_; i < alpha_start_; ++i) {
    vars_lowerbound[i] = min_acc_y_;
    vars_upperbound[i] = max_acc_y_;
  }

  for (auto i = alpha_start_; i < n_vars; ++i) {
    vars_lowerbound[i] = min_ang_acc_;
    vars_upperbound[i] = max_ang_acc_;
  }

  // setup constraints 
  auto n_constraints = mpc_step_ * state_dim;

  // constraints boundary conditions
  Dvector constraints_lowerbound(n_constraints), constraints_upperbound(n_constraints);

  for (auto i = 0; i < n_constraints; ++i) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start_] = x;
  constraints_lowerbound[y_start_] = y;
  constraints_lowerbound[psi_start_] = psi;
  constraints_lowerbound[vx_start_] = vx;
  constraints_lowerbound[vy_start_] = vy;
  constraints_lowerbound[cte_x_start_] = cte_x;
  constraints_lowerbound[cte_y_start_] = cte_y;

  constraints_upperbound[x_start_] = x;
  constraints_upperbound[y_start_] = y;
  constraints_upperbound[psi_start_] = psi;
  constraints_upperbound[vx_start_] = vx;
  constraints_upperbound[vy_start_] = vy;
  constraints_upperbound[cte_x_start_] = cte_x;
  constraints_upperbound[cte_y_start_] = cte_y;

  FG_eval fg_eval(trajs);
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

  return {
    answer[x_start_ + 1],  answer[y_start_ + 1], answer[psi_start_ + 1], 
    answer[vx_start_ + 1], answer[vy_start_ + 1], answer[w_start_ + 1],
    answer[cte_x_start_ + 1], answer[cte_y_start_ + 1],
    answer[ax_start_], answer[ay_start_], 
    answer[alpha_start_]
  };
}