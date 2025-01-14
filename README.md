# CarND Controls Quizzes

Quizzes for *Vehicle Models* and *Model Predictive Control* sections.

1. [Global Kinematic Model Quiz](./global_kinematic_model) - Implement the *Global Kinematic Model*.
2. [Polynomial Fitting Quiz](./polyfit) - Fit and evaluate polynomials.
3. [Mind The Line Quiz](./mpc_to_line) - Implement MPC and minimize cross track and orientation errors to a straight line trajectory.  See [this document](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/install_Ipopt_CppAD.md) for setup tips for executing the plotting code in the ```MPC.cpp``` solution file.

To do a quiz:

1. Go to quiz directory.
2. Make a build directory with `mkdir build`.
3. Change into the build directory, `cd build`.
4. Compile the project, `cmake .. && make`.

A solution for each quiz is presented in the solution directory.

## Dependencies

The *Global Kinematic Quiz* and *Polynomial Fitting* quizzes have all the dependencies in repo. For the *MPC* quiz
you'll have to install Ipopt and CppAD.  Please refer to [this document](https://github.com/udacity/CarND-MPC-Quizzes/blob/master/install_Ipopt_CppAD.md) for installation instructions.

```cpp
sudo apt-get install cppad
sudo apt-get install coinor-clp coinor-libclp-dev
sudo apt-get install coinor-libipopt-dev

# 22.04
sudo apt-get install python3.10-dev
# 20.04
sudo apt-get install python3.8-dev
```

```python
sudo apt-get install libeigen3-dev
sudo ln -s /usr/include/eigen3/Eigen /usr/include/Eigen
```