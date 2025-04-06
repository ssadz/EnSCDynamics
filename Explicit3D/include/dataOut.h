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

#ifndef DATAOUT_H
#define DATAOUT_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <chrono>
#include <iomanip>
#include "types.h"
#include "ElementHexN8.h"
#include "FSI_share_data.h"
namespace EnSC
{
	namespace DataOutType
	{
		enum DataType { scalar, vector, symmetric_tensor };
	}

	class DataOut {
	public:
		DataOut(const Types::VerticesAll<3>& pVerticesAll, const std::vector<Element_HexN8>& pTriaHex, Types::Real time);
		~DataOut() {}

		// ��ӽڵ�����
		void add_node_data(DataOutType::DataType dataType, const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>* vec, std::string name);

		// ��ӵ�Ԫ���ݣ������� dataType ����
		void add_ele_data(DataOutType::DataType dataType, const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>* vec, std::string name);
		void write_vtu(const std::string& fileName);
		void write_vtu_virtualParticles(const unsigned int& result_number, const FSI_share_data& fsi_data) const;
	private:
		const Types::VerticesAll<3>& opVerticesAll;
		const std::vector<Element_HexN8>& opTriaHex;
		Types::Real time;

		// �ڵ�����
		std::vector<DataOutType::DataType> opDataType;  // �ڵ���������
		std::vector<const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>*> opVec;
		std::vector<std::string> opName;

		// ��Ԫ����
		std::vector<DataOutType::DataType> opDataType_Ele;  // ��Ԫ��������
		std::vector<const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>*> opVec_Ele;
		std::vector<std::string> opName_Ele;
	};

}
#endif // DATAOUT_H
