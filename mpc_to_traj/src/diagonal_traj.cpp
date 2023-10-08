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

tuple<vector<double>, vector<double>> GetTrajPointsDiagonal(
    const int &start,
    const int &point_num,
    const double &distance,
    const int &sampling_num
){
    vector<double> x(point_num);
    vector<double> y(point_num);

    auto dx = distance / point_num;

    for(int i = 0; i < point_num; ++i)
    {
        x[i] = dx * i + start;
        y[i] = dx * i + start;
    }

    return std::make_tuple(x, y);
}

int main()
{
    const int point_num = 100;

    auto traj_points = GetTrajPointsDiagonal(0, point_num, 4, point_num);
    auto traj_x = std::get<0>(traj_points);
    auto traj_y = std::get<1>(traj_points);

    plt::plot(traj_x, traj_y);
    plt::grid(true); //show grid
    plt::show(); // show figure
}