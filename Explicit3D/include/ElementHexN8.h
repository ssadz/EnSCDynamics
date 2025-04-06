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

#ifndef ELEMENTHEXN8_H
#define ELEMENTHEXN8_H
#include "ElementBase.h"
namespace EnSC
{
	class Element_HexN8 : public ElementBase<3>
	{
	public:
		Element_HexN8();
		virtual ~Element_HexN8() {}
		Types::Real get_shapeFunctionValue(int pINode, const Types::Point<3>& pUnitPoint);
		Types::Real get_shapeFunctionDerivativeValue(int pINode, int pIDirection, const Types::Point<3>& pUnitPoint);
		std::array<Types::Real, 3> get_shapeFunctionDerivatives(int pINode, const std::array<Types::Real, 3>& pUnitPoint);
		static std::array<std::array<unsigned int, 4>, 6> Ele_face_indices;
	};

}
#endif // ELEMENTHEXN8_H
