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

#ifndef ELEMENTBASE_H
#define ELEMENTBASE_H
#include "types.h"
namespace EnSC
{
	template<int dim>
	class ElementBase
	{
	public:
		ElementBase(int nVertices);
		~ElementBase() {}
		const std::vector<Types::Vertex_index>& get_verticesIndex()const { return verticesIndex; }
		std::vector<Types::Vertex_index>& get_verticesIndex_writable() { return verticesIndex; }
		void set_verticesIndex(const std::vector<Types::Vertex_index>& pVerticesIndex);
		int get_PID() { return PID; }
		int get_MID() { return MID; }
		void set_PID(int pPID) { PID = pPID; }
		void set_MID(int pMID) { MID = pMID; }
		unsigned int get_nVertices() { return verticesIndex.size(); }
		virtual Types::Real get_shapeFunctionValue(int pINode, const Types::Point<dim>& pUnitPoint) = 0;
		virtual Types::Real get_shapeFunctionDerivativeValue(int pINode, int pIDirection, const Types::Point<dim>& pUnitPoint) = 0;
	protected:
		int PID;
		int MID;
		std::vector<Types::Vertex_index> verticesIndex;
	};

	//realize
	template<int dim>
	ElementBase<dim>::ElementBase(int nVertices)
	{
		verticesIndex.resize(nVertices);
	}

	template<int dim>
	void ElementBase<dim>::set_verticesIndex(const std::vector<Types::Vertex_index>& pVerticesIndex)
	{
		verticesIndex = pVerticesIndex;
	}

}
#endif // ELEMENTBASE_H
