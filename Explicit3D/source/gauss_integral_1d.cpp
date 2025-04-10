//----------------------------------------------------------------------------
//
// EnSC: 工程科学计算。EnSC的目的是提供
// 先进的工程科学计算数值工具包，不
// 限于有限元方法
//
// 基于Eigen3和STL开发。一些设计理念来自deal.II
//
// 这个文件是EnSC的一部分
//
// 作者: 盛文海
//
//----------------------------------------------------------------------------

#include "../include/gauss_integral_1d.h"

namespace EnSC
{
/**
 * @brief 构造函数实现
 * 
 * 根据所需的积分点数量初始化高斯积分点和权重
 * 支持1到6阶高斯积分，积分区间为[-1,1]
 * 
 * @param p_n_integral_points 积分点数量(1-6)
 */
Gauss_integral_1d::Gauss_integral_1d(int p_n_integral_points)
    :n_integral_points(p_n_integral_points)
{
    // 分配积分点和权重数组空间
    x.resize(n_integral_points);
    w.resize(n_integral_points);
    
    // 根据积分阶数设置相应的积分点和权重
    switch (n_integral_points)
    {
    case 1:  // 1点高斯积分（精确积分1阶多项式）
        x[0]=0.0; w[0]=2.0;
        break;

    case 2:  // 2点高斯积分（精确积分3阶多项式）
        x[0]=-0.577350269189626; w[0]=1.0;  // -1/√3
        x[1]=0.577350269189626;  w[1]=1.0;  // 1/√3
        break;

    case 3:  // 3点高斯积分（精确积分5阶多项式）
        x[0]=-0.774596669241483; w[0]=0.555555555555556;  // -√(3/5)
        x[1]=0.0;                w[1]=0.888888888888889;  // 0
        x[2]=0.774596669241483;  w[2]=0.555555555555556;  // √(3/5)
        break;

    case 4:  // 4点高斯积分（精确积分7阶多项式）
        x[0]=-0.861136311594053; w[0]=0.347854845137454;
        x[1]=-0.339981043584856; w[1]=0.652145154862546;
        x[2]=0.339981043584856;  w[2]=0.652145154862546;
        x[3]=0.861136311594053;  w[3]=0.347854845137454;
        break;

    case 5:  // 5点高斯积分（精确积分9阶多项式）
        x[0]=-0.906179845938664; w[0]=0.236926885056189;
        x[1]=-0.538469310105683; w[1]=0.478628670499366;
        x[2]=0.0;                w[2]=0.568888888888889;
        x[3]=0.538469310105683;  w[3]=0.478628670499366;
        x[4]=0.906179845938664;  w[4]=0.236926885056189;
        break;

    case 6:  // 6点高斯积分（精确积分11阶多项式）
        x[0]=-0.932469514203152; w[0]=0.171324492379170;
        x[1]=-0.661209386466265; w[1]=0.360761573048139;
        x[2]=-0.238619186083197; w[2]=0.467913934572691;
        x[3]=0.238619186083197;  w[3]=0.467913934572691;
        x[4]=0.661209386466265;  w[4]=0.360761573048139;
        x[5]=0.932469514203152;  w[5]=0.171324492379170;
        break;

    default:
        // 不支持的积分点数量，不进行初始化
        // 在实际应用中应该加入错误处理
        break;
    }
}

}
