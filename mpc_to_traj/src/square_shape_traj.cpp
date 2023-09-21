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

tuple<vector<double>, vector<double>> GetTrajPointsSquareShape(
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
        if(i < point_num / 4)
        {
            x[i] = dx * i + start;
            y[i] = 0;
        }
        else if(i < point_num / 2)
        {
            x[i] = dx * point_num / 4 + start;
            y[i] = 0 - dx * (i - point_num / 4);
        }
        else if(i < point_num * 3 / 4)
        {
            x[i] = dx * point_num / 4 + start - dx * (i - point_num / 2);
            y[i] = 0 - dx * point_num / 4;
        }
        else
        {
            x[i] = dx * point_num / 4 + start - dx * (point_num / 4);
            y[i] = 0 - dx * point_num / 4 + dx * (i - point_num * 3 / 4);
        }
    }

    return std::make_tuple(x, y);
}

int main()
{
    const int point_num = 100;

    auto traj_points = GetTrajPointsSquareShape(0, point_num, 40.0, point_num);
    auto traj_x = std::get<0>(traj_points);
    auto traj_y = std::get<1>(traj_points);

    plt::plot(traj_x, traj_y);
    plt::grid(true); //show grid
    plt::show(); // show figure
}