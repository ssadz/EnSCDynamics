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

#ifndef ELEMENTBASE_H
#define ELEMENTBASE_H
#include "types.h"
namespace EnSC
{
	/**
	 * @brief 单元基类模板
	 * 
	 * 该类是所有有限元单元的基类，定义了共同的属性和接口
	 * 使用模板参数dim表示单元的维度
	 * 
	 * @tparam dim 单元维度
	 */
	template<int dim>
	class ElementBase
	{
	public:
		/**
		 * @brief 构造函数
		 * @param nVertices 单元中的顶点数量
		 */
		ElementBase(int nVertices);
		
		/**
		 * @brief 析构函数
		 */
		~ElementBase() {}
		
		/**
		 * @brief 获取单元顶点索引的常量引用
		 * @return 顶点索引数组的常量引用
		 */
		const std::vector<Types::Vertex_index>& get_verticesIndex()const { return verticesIndex; }
		
		/**
		 * @brief 获取单元顶点索引的可写引用
		 * @return 顶点索引数组的可写引用
		 */
		std::vector<Types::Vertex_index>& get_verticesIndex_writable() { return verticesIndex; }
		
		/**
		 * @brief 设置单元的顶点索引
		 * @param pVerticesIndex 顶点索引数组
		 */
		void set_verticesIndex(const std::vector<Types::Vertex_index>& pVerticesIndex);
		
		/**
		 * @brief 获取单元的物理ID
		 * @return 物理ID
		 */
		int get_PID() { return PID; }
		
		/**
		 * @brief 获取单元的材料ID
		 * @return 材料ID
		 */
		int get_MID() { return MID; }
		
		/**
		 * @brief 设置单元的物理ID
		 * @param pPID 物理ID
		 */
		void set_PID(int pPID) { PID = pPID; }
		
		/**
		 * @brief 设置单元的材料ID
		 * @param pMID 材料ID
		 */
		void set_MID(int pMID) { MID = pMID; }
		
		/**
		 * @brief 获取单元的顶点数量
		 * @return 顶点数量
		 */
		unsigned int get_nVertices() { return verticesIndex.size(); }
		
		/**
		 * @brief 获取形函数值（纯虚函数）
		 * @param pINode 节点索引
		 * @param pUnitPoint 单位坐标系中的点坐标
		 * @return 形函数在该点的值
		 */
		virtual Types::Real get_shapeFunctionValue(int pINode, const Types::Point<dim>& pUnitPoint) = 0;
		
		/**
		 * @brief 获取形函数对指定方向的导数值（纯虚函数）
		 * @param pINode 节点索引
		 * @param pIDirection 方向索引
		 * @param pUnitPoint 单位坐标系中的点坐标
		 * @return 形函数对指定方向的导数值
		 */
		virtual Types::Real get_shapeFunctionDerivativeValue(int pINode, int pIDirection, const Types::Point<dim>& pUnitPoint) = 0;
	protected:
		int PID;  // 物理ID
		int MID;  // 材料ID
		std::vector<Types::Vertex_index> verticesIndex;  // 顶点索引数组
	};

	// 实现部分
	/**
	 * @brief ElementBase类构造函数实现
	 * @param nVertices 单元中的顶点数量
	 */
	template<int dim>
	ElementBase<dim>::ElementBase(int nVertices)
	{
		verticesIndex.resize(nVertices);
	}

	/**
	 * @brief 设置单元顶点索引的实现
	 * @param pVerticesIndex 顶点索引数组
	 */
	template<int dim>
	void ElementBase<dim>::set_verticesIndex(const std::vector<Types::Vertex_index>& pVerticesIndex)
	{
		verticesIndex = pVerticesIndex;
	}

}
#endif // ELEMENTBASE_H
