#ifndef MPC_H
#define MPC_H

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Core>
// #include "Eigen-3.3/Eigen/Core"

class MPC {
private:
  int mpc_step_;
  int x_start_, y_start_, psi_start_;
  int v_start_, w_start_; 
  int cte_start_, epsi_start_;
  int a_start_, alpha_start_;
  
  std::map<std::string, double> params_;

public:
  MPC();
  virtual ~MPC();

  void LoadParams(const std::map<std::string, double> &params);

  // Solve the model given an initial state.
  // Return the next state and actuations as a vector.
  std::vector<double> Solve(const Eigen::VectorXd &x0, 
                            const Eigen::VectorXd &coeffs);
};

#endif  // MPC_H
