//----------------------------------------------------------------------------
//
// EnSC: Engineering Scientific Computation. The purpose of EnSC is to provide
// advanced numerical toolkit for Engineering scientific computation, not
// limted to finit element method
//
// Based on Eigen3 , STL. Some design ideas come from deal.II
//
// This file is part of EnSC
//
// Authors: Sheng Wenhai
//
//----------------------------------------------------------------------------

#ifndef GAUSS_INTEGRAL_1D_H
#define GAUSS_INTEGRAL_1D_H
#include <vector>
namespace EnSC
{
class Gauss_integral_1d
{
public:
    Gauss_integral_1d(int p_n_integral_points);
    int get_n_integral_points()const {return n_integral_points;}
    double get_x(int p_i_integral_points)const {return x[p_i_integral_points];}
    double get_w(int p_i_integral_points)const {return w[p_i_integral_points];}

private:
    int n_integral_points;
    std::vector<double> x;
    std::vector<double> w;
};
}
#endif // GAUSS_INTEGRAL_1D_H
