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

#ifndef FSI_SHARE_DATA_H
#define FSI_SHARE_DATA_H
#include <vector>
#include <array>
#include"types.h"
struct FSI_share_data {
	std::vector<std::array<EnSC::Types::Real, 3>> FSI_virtualParticles_coordinates;
	std::vector<std::array<EnSC::Types::Real, 3>> FSI_virtualParticles_nodeForce;
	std::vector<std::array<EnSC::Types::Real, 3>> FSI_virtualParticles_velocity;
};

#endif // FSI_SHARE_DATA_H
