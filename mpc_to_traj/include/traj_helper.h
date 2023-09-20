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
    vector<double> x(point_num);
    vector<double> y(point_num);

    for(int i = 0; i < point_num; ++i)
    {
        x[i] = i + start;
        y[i] = 2;
    }

    return std::make_tuple(x, y);
}

/**
 * @brief Get the Traj Points Line With Option object
 * 
 * @param start_x : start point x 
 * @param total_point_num : total traj number
 * @param distance : distance from start point
 * @param sampling_num : get N points from full trajectory
 * @return tuple<vector<double>, vector<double>> 
 */
tuple<vector<double>, vector<double>> GetTrajPointsLineWithOption(
    const int &start_x,
    const int &total_point_num,
    const double &distance,
    const int &sampling_num
){
    vector<double> x(sampling_num);
    vector<double> y(sampling_num);

    auto dx = distance / total_point_num;

    for(int i = 0; i < sampling_num; ++i)
    {
        x[i] = dx * i + start_x;
        y[i] = 2;
    }

    return std::make_tuple(x, y);
}