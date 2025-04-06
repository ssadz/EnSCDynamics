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

#ifndef FSI_SHARE_DATA_H
#define FSI_SHARE_DATA_H
#include <vector>
#include <array>
#include"types.h"

/**
 * @brief 流固耦合共享数据结构
 * 
 * 该结构体用于存储流固耦合（Fluid-Structure Interaction, FSI）
 * 计算中虚拟粒子的坐标、节点力和速度数据，用于流体和结构之间的数据交换
 */
struct FSI_share_data {
	/// 虚拟粒子坐标数组，每个粒子有x,y,z三个坐标分量
	std::vector<std::array<EnSC::Types::Real, 3>> FSI_virtualParticles_coordinates;
	
	/// 作用在虚拟粒子上的节点力数组，每个力有x,y,z三个分量
	std::vector<std::array<EnSC::Types::Real, 3>> FSI_virtualParticles_nodeForce;
	
	/// 虚拟粒子速度数组，每个速度有x,y,z三个分量
	std::vector<std::array<EnSC::Types::Real, 3>> FSI_virtualParticles_velocity;
};

#endif // FSI_SHARE_DATA_H
