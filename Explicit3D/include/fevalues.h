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

#ifndef FEVALUES_H
#define FEVALUES_H
#include "types.h"
#include "gauss_integral_1d.h"
#include "ElementHexN8.h"

namespace EnSC
{
	/**
	 * @brief 有限元值类
	 * 
	 * 该类用于有限元计算中的形函数值和导数计算，包括雅可比矩阵计算
	 * 从单元参考空间到实际空间的映射、积分点和权重等
	 * 
	 * @tparam n_vertices_perEle 单元中顶点数量
	 * @tparam dim 空间维度
	 */
	template<int n_vertices_perEle, int dim>
	class FEValues
	{
	public:
		/**
		 * @brief 构造函数
		 * @param p_gauss_integral_1d 高斯积分对象
		 * @param p_verticesAll 所有节点坐标数组
		 */
		FEValues(const Gauss_integral_1d& p_gauss_integral_1d, const Types::VerticesAll<dim>& p_verticesAll);
		~FEValues() {}
		
		/**
		 * @brief 初始化单元类型相关数据
		 * @param p_eleType 单元类型枚举
		 */
		void init(Types::EleType p_eleType);
		
		/**
		 * @brief 重新初始化单元特定数据
		 * @param pEle 单元指针
		 */
		void reinit(const ElementBase<dim>* pEle);
		
		/**
		 * @brief 获取积分点处的形函数值
		 * @param iIntergralPoint 积分点索引
		 * @param iShapeFunc 形函数索引
		 * @return 形函数值
		 */
		Types::Real get_shape_value(int iIntergralPoint, Types::Vertex_index iShapeFunc) { return shape_func_vector_unit[iIntergralPoint](iShapeFunc); }
		
		/**
		 * @brief 获取积分点处参考单元上形函数导数值
		 * @param iIntergralPoint 积分点索引
		 * @param iShapeFunc 形函数索引
		 * @return 参考单元上形函数导数向量
		 */
		const Eigen::Matrix<Types::Real, dim, 1>& get_shape_derivatives_value_unit_node(int iIntergralPoint, Types::Vertex_index iShapeFunc) { return shape_derivatives_matrix_unit[iIntergralPoint].col(iShapeFunc); }
		
		/**
		 * @brief 获取积分点处实际单元上形函数导数值
		 * @param iIntergralPoint 积分点索引
		 * @param iShapeFunc 形函数索引
		 * @return 实际单元上形函数导数向量
		 */
		const Eigen::Matrix<Types::Real, dim, 1>& get_shape_derivatives_value_real_node(int iIntergralPoint, Types::Vertex_index iShapeFunc) { return shape_derivatives_matrix_real[iIntergralPoint].col(iShapeFunc); }
		
		/**
		 * @brief 获取积分点处实际单元上所有形函数导数矩阵
		 * @param iIntergralPoint 积分点索引
		 * @return 实际单元上形函数导数矩阵
		 */
		const Eigen::Matrix<Types::Real, dim, n_vertices_perEle>& get_shape_derivatives_value_real_matrix(int iIntergralPoint) { return shape_derivatives_matrix_real[iIntergralPoint]; }
		
		/**
		 * @brief 获取积分点处雅可比矩阵
		 * @param iIntergralPoint 积分点索引
		 * @return 雅可比矩阵
		 */
		const Eigen::Matrix<Types::Real, dim, dim>& get_jacobi(int iIntergralPoint) { return jacobi[iIntergralPoint]; }
		
		/**
		 * @brief 获取积分点处雅可比矩阵行列式
		 * @param iIntergralPoint 积分点索引
		 * @return 雅可比行列式值
		 */
		Types::Real get_jacobi_determinant(int iIntergralPoint) { return jacobi_determinant[iIntergralPoint]; }
		
		/**
		 * @brief 获取积分点处雅可比行列式与积分权重的乘积
		 * @param iIntergralPoint 积分点索引
		 * @return 雅可比行列式与积分权重的乘积
		 */
		Types::Real get_JxW(int iIntergralPoint) { return jacobi_determinant[iIntergralPoint] * integralW[iIntergralPoint]; }
		
		/**
		 * @brief 获取积分点处雅可比矩阵的逆
		 * @param iIntergralPoint 积分点索引
		 * @return 雅可比矩阵的逆
		 */
		const Eigen::Matrix<Types::Real, dim, dim>& get_inverse_jacobi(int iIntergralPoint) { return inverse_jacobi[iIntergralPoint]; }
		
		/**
		 * @brief 获取当前单元的节点坐标矩阵
		 * @return 节点坐标矩阵
		 */
		const Eigen::Matrix<Types::Real, n_vertices_perEle, dim>& get_xyzMatirx() { return xyz_matrix; }


	private:
		std::vector<Types::Point<dim>> integralPoints;        ///< 积分点在参考单元上的坐标
		std::vector<Types::Real> integralW;                   ///< 积分点权重
		std::vector<Types::Real> jacobi_determinant;          ///< 各积分点处雅可比矩阵行列式
		const Gauss_integral_1d& gauss_integral_1d;           ///< 一维高斯积分类引用
		const Types::VerticesAll<dim>& verticesAll;           ///< 所有节点坐标数组引用
		int n_integral_points;                                ///< 积分点总数
		int n_integral_points_1d;                             ///< 每个方向上的积分点数
		const ElementBase<dim>* ele = nullptr;                ///< 当前单元指针
		std::vector<Eigen::Matrix<Types::Real, 1, n_vertices_perEle>> shape_func_vector_unit;  ///< 参考单元上形函数值，每个积分点一个矩阵，1 x n_vertices_perEle
		std::vector<Eigen::Matrix<Types::Real, dim, n_vertices_perEle>> shape_derivatives_matrix_unit; ///< 参考单元上形函数导数，每个积分点一个矩阵，dim x n_vertices_perEle
		std::vector<Eigen::Matrix<Types::Real, dim, n_vertices_perEle>> shape_derivatives_matrix_real; ///< 实际单元上形函数导数，每个积分点一个矩阵，dim x n_vertices_perEle
		std::vector<Eigen::Matrix<Types::Real, dim, dim>> jacobi;                              ///< 各积分点处雅可比矩阵
		std::vector<Eigen::Matrix<Types::Real, dim, dim>> inverse_jacobi;                      ///< 各积分点处雅可比矩阵的逆
		Eigen::Matrix<Types::Real, n_vertices_perEle, dim> xyz_matrix;                         ///< 当前单元节点坐标矩阵，n_vertices_perEle x dim
	};

	//*********************************** 实现部分 ***********************************

	/**
	 * @brief 构造函数实现
	 */
	template<int n_vertices_perEle, int dim>
	FEValues<n_vertices_perEle, dim>::FEValues(const Gauss_integral_1d& p_gauss_integral_1d, const Types::VerticesAll<dim>& p_verticesAll)
		: gauss_integral_1d(p_gauss_integral_1d)
		, verticesAll(p_verticesAll)
	{
		// 根据维度确定积分点总数
		n_integral_points_1d = gauss_integral_1d.get_n_integral_points();
		if (dim == 3)
			n_integral_points = n_integral_points_1d * n_integral_points_1d * n_integral_points_1d;
		else if (dim == 2)
			n_integral_points = n_integral_points_1d * n_integral_points_1d;

		// 初始化各种容器
		integralPoints.resize(n_integral_points);
		integralW.resize(n_integral_points);
		jacobi_determinant.resize(n_integral_points);
		shape_func_vector_unit.resize(n_integral_points);
		shape_derivatives_matrix_unit.resize(n_integral_points);
		shape_derivatives_matrix_real.resize(n_integral_points);
		jacobi.resize(n_integral_points);
		inverse_jacobi.resize(n_integral_points);
	}

	/**
	 * @brief 初始化单元类型相关数据
	 * 
	 * 根据单元类型计算积分点坐标、权重、形函数值和导数
	 */
	template<int n_vertices_perEle, int dim>
	void FEValues<n_vertices_perEle, dim>::init(Types::EleType p_eleType)
	{
		Element_HexN8 tmp_HexN8;
		switch (p_eleType)
		{
		case Types::tetrahedron:
			break;

		case Types::hexN8:
			// 计算积分点坐标和权重
			for (int iIntegP = 0; iIntegP < n_integral_points; ++iIntegP) {
				auto& integralPoint = integralPoints[iIntegP];
				auto& integW = integralW[iIntegP];
				integW = 1.0;
				for (int i = 0; i < n_integral_points_1d; ++i)
					for (int j = 0; j < n_integral_points_1d; ++j)
						for (int k = 0; k < n_integral_points_1d; ++k)
						{
							integralPoint[0] = gauss_integral_1d.get_x(i);
							integW *= gauss_integral_1d.get_w(i);
							integralPoint[1] = gauss_integral_1d.get_x(j);
							integW *= gauss_integral_1d.get_w(j);
							integralPoint[2] = gauss_integral_1d.get_x(k);
							integW *= gauss_integral_1d.get_w(k);
						}
			}

			// 计算参考单元上的形函数值和导数
			for (int iIntegP = 0; iIntegP < n_integral_points; ++iIntegP)
			{
				const auto& integral_point = integralPoints[iIntegP];
				auto& iInteP_shape_func_uint = shape_func_vector_unit[iIntegP];
				auto& iInteP_shape_derivative_unit = shape_derivatives_matrix_unit[iIntegP];
				for (int iNode = 0; iNode < n_vertices_perEle; ++iNode)
				{
					iInteP_shape_func_uint(iNode) = tmp_HexN8.get_shapeFunctionValue(iNode, integral_point);
					for (int j = 0; j < dim; ++j)
						iInteP_shape_derivative_unit(j, iNode) = tmp_HexN8.get_shapeFunctionDerivativeValue(iNode, j, integral_point);
				}
			}
			break;

		case Types::wedge:
			break;

		case Types::pyramid:
			break;
		default:
			break;
		}
	}


	/**
	 * @brief 重新初始化单元特定数据
	 * 
	 * 根据当前单元计算节点坐标矩阵、雅可比矩阵及其逆、行列式，
	 * 以及实际单元上的形函数导数
	 */
	template<int n_vertices_perEle, int dim>
	void FEValues<n_vertices_perEle, dim>::reinit(const ElementBase<dim>* pEle)
	{
		// 保存单元指针
		ele = pEle;
		
		// 构建节点坐标矩阵
		const auto& tmpVertices = ele->get_verticesIndex();
		for (int i = 0; i < n_vertices_perEle; ++i)
		{
			const Types::Point<dim>& tmpPoint = verticesAll[tmpVertices[i]];
			xyz_matrix(i, 0) = tmpPoint[0];
			xyz_matrix(i, 1) = tmpPoint[1];
			xyz_matrix(i, 2) = tmpPoint[2];
		}
		
		// 对每个积分点计算雅可比矩阵及其逆
		for (int i = 0; i < n_integral_points; ++i)
		{
			// 计算雅可比矩阵: J = ∂x/∂ξ = ∑(∂Ni/∂ξ * xi)
			jacobi[i] = shape_derivatives_matrix_unit[i] * xyz_matrix;
			
			// 计算雅可比矩阵的逆和行列式
			inverse_jacobi[i] = jacobi[i].inverse();
			jacobi_determinant[i] = jacobi[i].determinant();
			
			// 计算实际单元上的形函数导数: ∂Ni/∂x = J^(-1) * ∂Ni/∂ξ
			auto& invJacobi = inverse_jacobi[i];
			for (int j = 0; j < n_vertices_perEle; ++j)
			{
				shape_derivatives_matrix_real[i].col(j) = invJacobi * shape_derivatives_matrix_unit[i].col(j);
			}
		}
	}

}
#endif // FEVALUES_H
