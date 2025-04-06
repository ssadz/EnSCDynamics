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

#include "../include/gauss_integral_1d.h"

namespace EnSC
{
Gauss_integral_1d::Gauss_integral_1d(int p_n_integral_points)
    :n_integral_points(p_n_integral_points)
{
    x.resize(n_integral_points);
    w.resize(n_integral_points);
    switch (n_integral_points)
    {
    case 1:
        x[0]=0.0; w[0]=2.0;
        break;

    case 2:
        x[0]=-0.577350269189626; w[0]=1.0;
        x[1]=0.577350269189626;  w[1]=1.0;
        break;

    case 3:
        x[0]=-0.774596669241483; w[0]=0.555555555555556;
        x[1]=0.0;                w[1]=0.888888888888889;
        x[2]=0.774596669241483;  w[2]=0.555555555555556;
        break;

    case 4:
        x[0]=-0.861136311594053; w[0]=0.347854845137454;
        x[1]=-0.339981043584856; w[1]=0.652145154862546;
        x[2]=0.339981043584856;  w[2]=0.652145154862546;
        x[3]=0.861136311594053;  w[3]=0.347854845137454;
        break;

    case 5:
        x[0]=-0.906179845938664; w[0]=0.236926885056189;
        x[1]=-0.538469310105683; w[1]=0.478628670499366;
        x[2]=0.0;                w[2]=0.568888888888889;
        x[3]=0.538469310105683;  w[3]=0.478628670499366;
        x[4]=0.906179845938664;  w[4]=0.236926885056189;
        break;

    case 6:
        x[0]=-0.932469514203152; w[0]=0.171324492379170;
        x[1]=-0.661209386466265; w[1]=0.360761573048139;
        x[2]=-0.238619186083197; w[2]=0.467913934572691;
        x[3]=0.238619186083197;  w[3]=0.467913934572691;
        x[4]=0.661209386466265;  w[4]=0.360761573048139;
        x[5]=0.932469514203152;  w[5]=0.171324492379170;
        break;

    default:
        break;
    }
}

}
