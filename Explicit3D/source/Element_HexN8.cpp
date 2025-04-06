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

#include "../include/ElementHexN8.h"

namespace EnSC {
	// 定义六面体单元六个面的节点索引
	// 每一行代表一个面，包含4个节点的局部索引
	// 面的编号顺序：下(z-)、上(z+)、前(y-)、后(y+)、左(x-)、右(x+)
	std::array<std::array<unsigned int, 4>, 6> Element_HexN8::Ele_face_indices = { {{3,0,4,7},  // 下面 z-
																		 {1,2,6,5},  // 上面 z+
																		 {0,1,5,4},  // 前面 y-
																		 {2,3,7,6},  // 后面 y+
																		 {1,0,3,2},  // 左面 x-
																		 {4,5,6,7}} };  // 右面 x+
	/**
	 * @brief 构造函数实现
	 * 
	 * 调用基类构造函数，传入顶点数量8
	 */
	Element_HexN8::Element_HexN8()
		:ElementBase(8) {
	}

	/**
	 * @brief 计算形函数值实现
	 * 
	 * 根据参考单元上的节点坐标和给定点坐标计算形函数值
	 * 使用标准三线性形函数：N_i = (1/8)*(1+ξ_i*ξ)*(1+η_i*η)*(1+ζ_i*ζ)
	 * 
	 * @param pINode 节点索引(0-7)
	 * @param pUnitPoint 参考单元中的点坐标
	 * @return 形函数在该点的值
	 */
	Types::Real Element_HexN8::get_shapeFunctionValue(int pINode, const Types::Point<3>& pUnitPoint) {
		// 定义参考单元上节点的坐标
		Types::Point<3> vertexCoor_unit;
		switch (pINode) {
		case 0:  // 节点0: (-1,-1,-1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 1:  // 节点1: (1,-1,-1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 2:  // 节点2: (1,1,-1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 3:  // 节点3: (-1,1,-1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 4:  // 节点4: (-1,-1,1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 5:  // 节点5: (1,-1,1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 6:  // 节点6: (1,1,1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 7:  // 节点7: (-1,1,1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		default:
			break;
		}
		// 计算三线性形函数值: N_i = (1/8)*(1+ξ_i*ξ)*(1+η_i*η)*(1+ζ_i*ζ)
		return 1.0 / 8.0 * (1 + vertexCoor_unit[0] * pUnitPoint[0]) * (1 + vertexCoor_unit[1] * pUnitPoint[1]) * (1 + vertexCoor_unit[2] * pUnitPoint[2]);
	}

	/**
	 * @brief 计算形函数导数值实现
	 * 
	 * 根据参考单元上的节点坐标和给定点坐标计算形函数对指定方向的导数值
	 * 
	 * @param pINode 节点索引(0-7)
	 * @param pIDirection 求导方向(0,1,2 分别对应 x,y,z)
	 * @param pUnitPoint 参考单元中的点坐标
	 * @return 形函数在该点对指定方向的导数值
	 */
	Types::Real Element_HexN8::get_shapeFunctionDerivativeValue(int pINode, int pIDirection, const Types::Point<3>& pUnitPoint) {
		// 定义参考单元上节点的坐标
		Types::Point<3> vertexCoor_unit;
		switch (pINode) {
		case 0:  // 节点0: (-1,-1,-1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 1:  // 节点1: (1,-1,-1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 2:  // 节点2: (1,1,-1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 3:  // 节点3: (-1,1,-1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 4:  // 节点4: (-1,-1,1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 5:  // 节点5: (1,-1,1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 6:  // 节点6: (1,1,1)
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 7:  // 节点7: (-1,1,1)
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		default:
			break;
		}
		
		// 根据方向计算形函数的偏导数
		if (pIDirection == 0)
			// 对x方向的偏导数: ∂N_i/∂ξ = (1/8)*ξ_i*(1+η_i*η)*(1+ζ_i*ζ)
			return 1.0 / 8.0 * vertexCoor_unit[0] * (1 + vertexCoor_unit[1] * pUnitPoint[1]) * (1 + vertexCoor_unit[2] * pUnitPoint[2]);
		else if (pIDirection == 1)
			// 对y方向的偏导数: ∂N_i/∂η = (1/8)*(1+ξ_i*ξ)*η_i*(1+ζ_i*ζ)
			return 1.0 / 8.0 * (1 + vertexCoor_unit[0] * pUnitPoint[0]) * vertexCoor_unit[1] * (1 + vertexCoor_unit[2] * pUnitPoint[2]);
		else
			// 对z方向的偏导数: ∂N_i/∂ζ = (1/8)*(1+ξ_i*ξ)*(1+η_i*η)*ζ_i
			return 1.0 / 8.0 * (1 + vertexCoor_unit[0] * pUnitPoint[0]) * (1 + vertexCoor_unit[1] * pUnitPoint[1]) * vertexCoor_unit[2];
	}

	/**
	 * @brief 计算形函数对所有方向的导数值
	 * 
	 * 同时计算形函数对x,y,z三个方向的导数，返回导数数组
	 * 
	 * @param pINode 节点索引(0-7)
	 * @param pUnitPoint 参考单元中的点坐标
	 * @return 包含形函数对x,y,z三个方向导数的数组
	 */
	std::array<Types::Real, 3> Element_HexN8::get_shapeFunctionDerivatives(int pINode, const std::array<double, 3>& pUnitPoint)
	{
		std::array<Types::Real, 3> derivatives;
		// 计算形函数对三个方向的导数
		for (int i = 0; i < 3; ++i) {
			derivatives[i] = get_shapeFunctionDerivativeValue(pINode, i, pUnitPoint);
		}
		return derivatives;
	}
}
