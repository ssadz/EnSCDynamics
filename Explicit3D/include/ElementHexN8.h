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

#ifndef ELEMENTHEXN8_H
#define ELEMENTHEXN8_H
#include "ElementBase.h"
namespace EnSC
{
	/**
	 * @brief 8节点六面体单元类
	 * 
	 * 该类实现了8节点六面体单元的形函数计算和形函数导数计算
	 * 继承自三维单元基类ElementBase<3>
	 */
	class Element_HexN8 : public ElementBase<3>
	{
	public:
		/**
		 * @brief 默认构造函数
		 */
		Element_HexN8();
		
		/**
		 * @brief 虚析构函数
		 */
		virtual ~Element_HexN8() {}
		
		/**
		 * @brief 获取形函数值
		 * @param pINode 节点索引(0-7)
		 * @param pUnitPoint 单位坐标系中的点坐标
		 * @return 形函数在该点的值
		 */
		Types::Real get_shapeFunctionValue(int pINode, const Types::Point<3>& pUnitPoint);
		
		/**
		 * @brief 获取形函数对单一方向的导数值
		 * @param pINode 节点索引(0-7)
		 * @param pIDirection 方向索引(0-2)，对应x,y,z方向
		 * @param pUnitPoint 单位坐标系中的点坐标
		 * @return 形函数对指定方向的导数值
		 */
		Types::Real get_shapeFunctionDerivativeValue(int pINode, int pIDirection, const Types::Point<3>& pUnitPoint);
		
		/**
		 * @brief 获取形函数对所有方向的导数值
		 * @param pINode 节点索引(0-7)
		 * @param pUnitPoint 单位坐标系中的点坐标
		 * @return 包含形函数对x,y,z三个方向导数的数组
		 */
		std::array<Types::Real, 3> get_shapeFunctionDerivatives(int pINode, const std::array<Types::Real, 3>& pUnitPoint);
		
		/**
		 * @brief 六面体单元的面节点索引定义
		 * 
		 * 静态数组，定义了六面体单元六个面上节点的局部索引
		 * 每个面由4个节点定义，共6个面
		 */
		static std::array<std::array<unsigned int, 4>, 6> Ele_face_indices;
	};

}
#endif // ELEMENTHEXN8_H
