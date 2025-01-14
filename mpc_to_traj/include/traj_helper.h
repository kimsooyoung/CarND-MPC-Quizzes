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

// L shape trajectory

/**
 * @brief Get the Traj Points L Shape object
 * 
 * @param start x start location
 * @param point_num total points
 * @param distance full path distance
 * @param sampling_num subsampling numbers
 * @return tuple<vector<double>, vector<double>> 
 */
tuple<vector<double>, vector<double>> GetTrajPointsLShape(
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
        if(i < point_num / 2)
        {
            x[i] = dx * i + start;
            y[i] = 10;
        }
        else
        {
            x[i] = dx * point_num / 2 + start;
            y[i] = 10 - dx * (i - point_num / 2);
        }
    }

    return std::make_tuple(x, y);
}

// square shape trajectory
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

// diagonal shape trajectory
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
