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
		using Real = double;
		using Vertex_index = unsigned int;
		using Dof_index = unsigned int;

		enum EleType
		{
			tetrahedron, hexN8, wedge, pyramid
		};

		template<int dim> using Point = std::array<Types::Real, dim>;

		template<int dim> using VerticesAll = std::vector<Point<dim>>;
	}

	constexpr Types::Real zero = (Types::Real)0.0, one = (Types::Real)1.0, two = (Types::Real)2.0,
		three = (Types::Real)3.0, four = (Types::Real)4.0, five = (Types::Real)5.0,
		six = (Types::Real)6.0, seven = (Types::Real)7.0, eight = (Types::Real)8.0,
		nine = (Types::Real)9.0, billion = (Types::Real)1.0e9;
	constexpr Types::Real qrt = one / four;
}

#endif // TYPES_H
