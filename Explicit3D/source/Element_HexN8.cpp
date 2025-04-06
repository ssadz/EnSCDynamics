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

#include "../include/ElementHexN8.h"

namespace EnSC {
	std::array<std::array<unsigned int, 4>, 6> Element_HexN8::Ele_face_indices = { {{3,0,4,7},
																		 {1,2,6,5},
																		 {0,1,5,4},
																		 {2,3,7,6},
																		 {1,0,3,2},
																		 {4,5,6,7}} };
	Element_HexN8::Element_HexN8()
		:ElementBase(8) {
	}

	Types::Real Element_HexN8::get_shapeFunctionValue(int pINode, const Types::Point<3>& pUnitPoint) {
		Types::Point<3> vertexCoor_unit;
		switch (pINode) {
		case 0:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 1:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 2:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 3:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 4:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 5:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 6:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 7:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		default:
			break;
		}
		return 1.0 / 8.0 * (1 + vertexCoor_unit[0] * pUnitPoint[0]) * (1 + vertexCoor_unit[1] * pUnitPoint[1]) * (1 + vertexCoor_unit[2] * pUnitPoint[2]);
	}

	Types::Real Element_HexN8::get_shapeFunctionDerivativeValue(int pINode, int pIDirection, const Types::Point<3>& pUnitPoint) {
		Types::Point<3> vertexCoor_unit;
		switch (pINode) {
		case 0:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 1:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 2:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 3:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = -1.0;
			break;
		case 4:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 5:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = -1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 6:
			vertexCoor_unit[0] = 1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		case 7:
			vertexCoor_unit[0] = -1.0;
			vertexCoor_unit[1] = 1.0;
			vertexCoor_unit[2] = 1.0;
			break;
		default:
			break;
		}
		if (pIDirection == 0)
			return 1.0 / 8.0 * vertexCoor_unit[0] * (1 + vertexCoor_unit[1] * pUnitPoint[1]) * (1 + vertexCoor_unit[2] * pUnitPoint[2]);
		else if (pIDirection == 1)
			return 1.0 / 8.0 * (1 + vertexCoor_unit[0] * pUnitPoint[0]) * vertexCoor_unit[1] * (1 + vertexCoor_unit[2] * pUnitPoint[2]);
		else
			return 1.0 / 8.0 * (1 + vertexCoor_unit[0] * pUnitPoint[0]) * (1 + vertexCoor_unit[1] * pUnitPoint[1]) * vertexCoor_unit[2];
	}

	std::array<Types::Real, 3> Element_HexN8::get_shapeFunctionDerivatives(int pINode, const std::array<double, 3>& pUnitPoint)
	{
		std::array<Types::Real, 3> derivatives;
		for (int i = 0; i < 3; ++i) {
			derivatives[i] = get_shapeFunctionDerivativeValue(pINode, i, pUnitPoint);
		}
		return derivatives;
	}
}
