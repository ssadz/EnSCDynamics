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

#ifndef TYPES_H
#define TYPES_H
#include <array>
#include <vector>
#include <Eigen/Dense>
#include <locale>
namespace EnSC
{

	namespace Types
	{
		// 使用双精度浮点数作为实数类型
		using Real = double;
		// 顶点索引类型定义
		using Vertex_index = unsigned int;
		// 自由度索引类型定义
		using Dof_index = unsigned int;

		// 单元类型枚举
		enum EleType
		{
			tetrahedron, // 四面体单元
			hexN8,       // 8节点六面体单元
			wedge,       // 楔形单元
			pyramid      // 金字塔单元
		};

		// 定义dim维空间中的点类型
		template<int dim> using Point = std::array<Types::Real, dim>;

		// 定义dim维空间中所有顶点的集合类型
		template<int dim> using VerticesAll = std::vector<Point<dim>>;
	}

	// 定义常用数值常量，使用constexpr提高计算效率
	constexpr Types::Real zero = (Types::Real)0.0, one = (Types::Real)1.0, two = (Types::Real)2.0,
		three = (Types::Real)3.0, four = (Types::Real)4.0, five = (Types::Real)5.0,
		six = (Types::Real)6.0, seven = (Types::Real)7.0, eight = (Types::Real)8.0,
		nine = (Types::Real)9.0, billion = (Types::Real)1.0e9;
	constexpr Types::Real qrt = one / four; // 四分之一常量，用于计算
}

#endif // TYPES_H
