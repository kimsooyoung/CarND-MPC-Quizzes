#pragma once

#include <vector>
#include <cmath>
#include <tuple>

using std::vector;
using std::tuple;
using std::sin;
using std::cos;

tuple<vector<double>, vector<double>> GetTrajPointsCirc(
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

tuple<vector<double>, vector<double>> GetTrajPointsLine(
    const int &start,
    const int &point_num,
    const int &sampling_num
){

    const double pi = M_PI;
    vector<double> x(point_num);
    vector<double> y(point_num);

    for(int i = 0; i < point_num; ++i)
    {
        x[i] = i + start;
        y[i] = 2;
    }

    return std::make_tuple(x, y);
}