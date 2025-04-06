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

#ifndef GAUSS_INTEGRAL_1D_H
#define GAUSS_INTEGRAL_1D_H
#include <vector>
namespace EnSC
{
/**
 * @brief 一维高斯积分类
 * 
 * 该类用于实现一维高斯积分，提供高斯积分点和权重，
 * 支持不同阶数的高斯积分
 */
class Gauss_integral_1d
{
public:
    /**
     * @brief 构造函数
     * @param p_n_integral_points 高斯积分点的数量（积分阶数）
     */
    Gauss_integral_1d(int p_n_integral_points);
    
    /**
     * @brief 获取高斯积分点的数量
     * @return 积分点数量
     */
    int get_n_integral_points()const {return n_integral_points;}
    
    /**
     * @brief 获取指定索引的高斯积分点坐标
     * @param p_i_integral_points 积分点索引
     * @return 积分点坐标（在[-1,1]区间内）
     */
    double get_x(int p_i_integral_points)const {return x[p_i_integral_points];}
    
    /**
     * @brief 获取指定索引的高斯积分点权重
     * @param p_i_integral_points 积分点索引
     * @return 积分点权重
     */
    double get_w(int p_i_integral_points)const {return w[p_i_integral_points];}

private:
    int n_integral_points;    ///< 高斯积分点的数量
    std::vector<double> x;    ///< 高斯积分点坐标数组
    std::vector<double> w;    ///< 高斯积分点权重数组
};
}
#endif // GAUSS_INTEGRAL_1D_H
