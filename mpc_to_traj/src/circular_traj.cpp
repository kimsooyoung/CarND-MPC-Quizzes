#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <Eigen/QR>
#include <matplotlibcpp.h>

#include "helpers.h"

using Eigen::VectorXd;
using std::vector;
using std::tuple;
using std::sin;
using std::cos;
using std::cout;

namespace plt = matplotlibcpp;

tuple<vector<double>, vector<double>> GetTrajPoints(
    const int &start,
    const int &point_num,
    const int &sampling_num
){

    const double pi = M_PI;
    vector<double> x(point_num);
    vector<double> y(point_num);

    for(int i = 0; i < point_num; ++i)
    {
        x[i] = cos(2 * pi * (i + start % sampling_num) / sampling_num - pi / 2);
        y[i] = sin(2 * pi * (i + start % sampling_num) / sampling_num - pi / 2 ) + 1;
    }

    return std::make_tuple(x, y);
}

int main()
{
    const double pi = M_PI;
    const int point_num = 5;
    const int iter_num = 200;
    
    vector<double> x(iter_num + 1);
    vector<double> y(iter_num + 1);

    for(int i = 0; i < iter_num + 1; ++i)
    {
        x[i] = cos(2 * pi * i / iter_num);
        y[i] = sin(2 * pi * i / iter_num) + 1;
    }

    // GetTrajPoints(start, point_num, sampling_num)
    auto traj_points = GetTrajPoints(0, point_num, iter_num);
    auto traj_x = std::get<0>(traj_points);
    auto traj_y = std::get<1>(traj_points);

    for (auto x : traj_x)
        std::cout << "traj_x : " << x << std::endl;
    for (auto y : traj_y)
        std::cout << "traj_y : " << y << std::endl;

    VectorXd x_veh(point_num);
    VectorXd y_veh(point_num);
    for(long unsigned int i = 0; i < point_num; ++i)
    {
        x_veh[i] = traj_x[i];
        y_veh[i] = traj_y[i];
    }

    auto coeffs = polyfit(x_veh, y_veh, 1);
    std::cout << coeffs << std::endl;

    plt::plot(x, y, "r--"); //plot the x,y
    plt::plot(traj_x, traj_y);
    plt::axis("scaled"); // equal axis scale
    plt::grid(true); //show grid
    plt::show(); // show figure
}