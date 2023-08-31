#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <matplotlibcpp.h>

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
    const int iter_num = 50;
    
    vector<double> x(iter_num + 1);
    vector<double> y(iter_num + 1);

    for(int i = 0; i < iter_num + 1; ++i)
    {
        x[i] = cos(2 * pi * i / iter_num);
        y[i] = sin(2 * pi * i / iter_num) + 1;
    }

    // GetTrajPoints(start, point_num, sampling_num)
    auto traj_points = GetTrajPoints(25, 6, iter_num);
    auto traj_x = std::get<0>(traj_points);
    auto traj_y = std::get<1>(traj_points);

    plt::plot(x, y, "r--"); //plot the x,y
    plt::plot(traj_x, traj_y);
    plt::grid(true); //show grid
    plt::show(); // show figure
}