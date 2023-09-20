#ifndef MPC_H
#define MPC_H

#include <map>
#include <tuple>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Core>
// #include "Eigen-3.3/Eigen/Core"

using std::vector;
using std::tuple;

class MPCQD {
private:
  int mpc_step_;
  int x_start_, y_start_, psi_start_;
  int vx_start_, vy_start_, w_start_; 
  int cte_x_start_, cte_y_start_;
  int ax_start_, ay_start_, alpha_start_;

  double min_acc_x_, max_acc_x_;
  double min_acc_y_, max_acc_y_;
  double min_ang_acc_, max_ang_acc_;
  
  std::map<std::string, double> params_;

public:
  MPCQD();
  virtual ~MPCQD();

  void LoadParams(const std::map<std::string, double> &params);

  // Solve the model given an initial state.
  // Return the next state and actuations as a vector.
  std::vector<double> Solve(const Eigen::VectorXd &x0, 
                            const std::tuple<vector<double>, vector<double>> &trajs);
};

#endif  // MPC_H
