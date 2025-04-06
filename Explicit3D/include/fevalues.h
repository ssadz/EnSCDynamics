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

#ifndef FEVALUES_H
#define FEVALUES_H
#include "types.h"
#include "gauss_integral_1d.h"
#include "ElementHexN8.h"

namespace EnSC
{
	template<int n_vertices_perEle, int dim>
	class FEValues
	{
	public:
		FEValues(const Gauss_integral_1d& p_gauss_integral_1d, const Types::VerticesAll<dim>& p_verticesAll);
		~FEValues() {}
		void init(Types::EleType p_eleType);
		void reinit(const ElementBase<dim>* pEle);
		Types::Real get_shape_value(int iIntergralPoint, Types::Vertex_index iShapeFunc) { return shape_func_vector_unit[iIntergralPoint](iShapeFunc); }
		const Eigen::Matrix<Types::Real, dim, 1>& get_shape_derivatives_value_unit_node(int iIntergralPoint, Types::Vertex_index iShapeFunc) { return shape_derivatives_matrix_unit[iIntergralPoint].col(iShapeFunc); }
		const Eigen::Matrix<Types::Real, dim, 1>& get_shape_derivatives_value_real_node(int iIntergralPoint, Types::Vertex_index iShapeFunc) { return shape_derivatives_matrix_real[iIntergralPoint].col(iShapeFunc); }
		const Eigen::Matrix<Types::Real, dim, n_vertices_perEle>& get_shape_derivatives_value_real_matrix(int iIntergralPoint) { return shape_derivatives_matrix_real[iIntergralPoint]; }
		const Eigen::Matrix<Types::Real, dim, dim>& get_jacobi(int iIntergralPoint) { return jacobi[iIntergralPoint]; }
		Types::Real get_jacobi_determinant(int iIntergralPoint) { return jacobi_determinant[iIntergralPoint]; }
		Types::Real get_JxW(int iIntergralPoint) { return jacobi_determinant[iIntergralPoint] * integralW[iIntergralPoint]; }
		const Eigen::Matrix<Types::Real, dim, dim>& get_inverse_jacobi(int iIntergralPoint) { return inverse_jacobi[iIntergralPoint]; }
		const Eigen::Matrix<Types::Real, n_vertices_perEle, dim>& get_xyzMatirx() { return xyz_matrix; }


	private:
		std::vector<Types::Point<dim>> integralPoints;
		std::vector<Types::Real> integralW;
		std::vector<Types::Real> jacobi_determinant;
		const Gauss_integral_1d& gauss_integral_1d;
		const Types::VerticesAll<dim>& verticesAll;
		int n_integral_points;
		int n_integral_points_1d;
		const ElementBase<dim>* ele = nullptr;
		std::vector<Eigen::Matrix<Types::Real, 1, n_vertices_perEle>> shape_func_vector_unit;//1 x n_vertices_perEle
		std::vector<Eigen::Matrix<Types::Real, dim, n_vertices_perEle>> shape_derivatives_matrix_unit;//3 x n_vertices_perEle
		std::vector<Eigen::Matrix<Types::Real, dim, n_vertices_perEle>> shape_derivatives_matrix_real;//3 x n_vertices_perEle
		std::vector<Eigen::Matrix<Types::Real, dim, dim>> jacobi;
		std::vector<Eigen::Matrix<Types::Real, dim, dim>> inverse_jacobi;
		Eigen::Matrix<Types::Real, n_vertices_perEle, dim> xyz_matrix;//n_vertices_perEle x 3
	};

	//***********************************Realize***********************************

	template<int n_vertices_perEle, int dim>
	FEValues<n_vertices_perEle, dim>::FEValues(const Gauss_integral_1d& p_gauss_integral_1d, const Types::VerticesAll<dim>& p_verticesAll)
		: gauss_integral_1d(p_gauss_integral_1d)
		, verticesAll(p_verticesAll)
	{
		n_integral_points_1d = gauss_integral_1d.get_n_integral_points();
		if (dim == 3)
			n_integral_points = n_integral_points_1d * n_integral_points_1d * n_integral_points_1d;
		else if (dim == 2)
			n_integral_points = n_integral_points_1d * n_integral_points_1d;

		integralPoints.resize(n_integral_points);
		integralW.resize(n_integral_points);
		jacobi_determinant.resize(n_integral_points);
		shape_func_vector_unit.resize(n_integral_points);
		shape_derivatives_matrix_unit.resize(n_integral_points);
		shape_derivatives_matrix_real.resize(n_integral_points);
		jacobi.resize(n_integral_points);
		inverse_jacobi.resize(n_integral_points);
	}

	template<int n_vertices_perEle, int dim>
	void FEValues<n_vertices_perEle, dim>::init(Types::EleType p_eleType)
	{
		//shape_func_vector_unit,shape_derivatives_matrix_unit
		Element_HexN8 tmp_HexN8;
		switch (p_eleType)
		{
		case Types::tetrahedron:
			break;

		case Types::hexN8:
			//dim integralPoints and integralW
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

			//shape_func_vector_unit,shape_derivatives_matrix_unit
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


	template<int n_vertices_perEle, int dim>
	void FEValues<n_vertices_perEle, dim>::reinit(const ElementBase<dim>* pEle)
	{
		//ele
		ele = pEle;
		//xyz_matrix
		const auto& tmpVertices = ele->get_verticesIndex();
		for (int i = 0; i < n_vertices_perEle; ++i)
		{
			const Types::Point<dim>& tmpPoint = verticesAll[tmpVertices[i]];
			xyz_matrix(i, 0) = tmpPoint[0];
			xyz_matrix(i, 1) = tmpPoint[1];
			xyz_matrix(i, 2) = tmpPoint[2];
		}
		for (int i = 0; i < n_integral_points; ++i)
		{
			//jacobi
			//inverse_jacobi
			jacobi[i] = shape_derivatives_matrix_unit[i] * xyz_matrix;
			inverse_jacobi[i] = jacobi[i].inverse();
			jacobi_determinant[i] = jacobi[i].determinant();
			//shape_derivatives_matrix_real
			auto& invJacobi = inverse_jacobi[i];
			for (int j = 0; j < n_vertices_perEle; ++j)
			{
				shape_derivatives_matrix_real[i].col(j) = invJacobi * shape_derivatives_matrix_unit[i].col(j);
			}
		}
	}

}
#endif // FEVALUES_H
