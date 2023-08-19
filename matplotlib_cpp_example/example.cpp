#include <iostream>
#include <vector>
#include <cmath>
#include <matplotlibcpp.h>

using std::vector;
using std::sin;
using std::cout;

namespace plt = matplotlibcpp;

int main()
{
    const double pi = M_PI;
    
    vector<double> x(21);
    vector<double> y(21);

    for(int i = 1; i < x.capacity(); ++i)
    {
        x[i] = x[i-1] + pi/10;
        y[i] = sin(x[i]);
    }

    plt::plot(x,y); //plot the x,y
    plt::grid(true); //show grid
    plt::show(); // show figure
}