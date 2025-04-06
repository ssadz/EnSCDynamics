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

#include "../include/dataOut.h"
#include <sys/stat.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <sys/stat.h>
#include <vtkDoubleArray.h>
#ifdef _WIN32
#include <direct.h> // <<<--- 添加这个头文件 for _mkdir
#endif
namespace EnSC
{
	DataOut::DataOut(const Types::VerticesAll<3>& pVerticesAll, const std::vector<Element_HexN8>& pTriaHex, Types::Real time)
		: opVerticesAll(pVerticesAll),
		opTriaHex(pTriaHex),
		time(time) {
	}


	void DataOut::add_node_data(DataOutType::DataType dataType, const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>* vec, std::string name) {
		opDataType.emplace_back(dataType);
		opVec.emplace_back(vec);
		opName.emplace_back(name);
	}

	void DataOut::add_ele_data(DataOutType::DataType dataType, const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>* vec, std::string name) {
		opDataType_Ele.emplace_back(dataType);
		opVec_Ele.emplace_back(vec);
		opName_Ele.emplace_back(name);
	}

	void DataOut::write_vtu(const std::string& fileName) {
		// Check and create output directory
		struct stat info;
		if (stat("ResultData", &info) != 0 || !(info.st_mode & S_IFDIR)) {
			system("mkdir ResultData");
		}

		// Construct full output file path
		std::string filepath = "ResultData/" + fileName;

		// Create a new unstructured grid
		auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

		// Create and set points
		auto points = vtkSmartPointer<vtkPoints>::New();
		for (const auto& vertex : opVerticesAll) {
			points->InsertNextPoint(vertex[0], vertex[1], vertex[2]); // Now using z-coordinate
		}
		grid->SetPoints(points);

		// Create and set cells as VTK_HEXAHEDRON
		for (const auto& ele : opTriaHex) {
			const auto& vertices = ele.get_verticesIndex();

			vtkIdType pts[8];
			for (vtkIdType i = 0; i < 8; ++i) {
				pts[i] = vertices[i];
			}

			// Set cell type to VTK_HEXAHEDRON
			grid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
		}

		// Add point data
		for (size_t i = 0; i < opDataType.size(); ++i) {
			vtkSmartPointer<vtkFloatArray> dataArray = vtkSmartPointer<vtkFloatArray>::New();
			dataArray->SetName(opName[i].c_str());
			if (opDataType[i] == DataOutType::scalar) {
				dataArray->SetNumberOfComponents(1); // Scalar
				dataArray->SetNumberOfTuples(opVerticesAll.size()); // Pre-allocate space
#pragma omp parallel for
				for (Types::Vertex_index j = 0; j < opVerticesAll.size(); ++j) {
					dataArray->SetValue(j, (*opVec[i])[j]); // Fill data
				}
			}
			else if (opDataType[i] == DataOutType::vector) {
				dataArray->SetNumberOfComponents(3); // Vector has 3 components
				dataArray->SetNumberOfTuples(opVerticesAll.size()); // Pre-allocate space
#pragma omp parallel for
				for (Types::Vertex_index j = 0; j < opVerticesAll.size(); ++j) {
					dataArray->SetTuple3(j, (*opVec[i])[3 * j], (*opVec[i])[3 * j + 1], (*opVec[i])[3 * j + 2]); // Fill data
				}
			}
			// If need to handle symmetric tensor data at nodes, can add here
			grid->GetPointData()->AddArray(dataArray);
		}

		// Add cell data
		for (size_t i = 0; i < opName_Ele.size(); ++i) {
			vtkSmartPointer<vtkFloatArray> dataArray = vtkSmartPointer<vtkFloatArray>::New();
			dataArray->SetName(opName_Ele[i].c_str());
			if (opDataType_Ele[i] == DataOutType::scalar) {
				// Scalar data
				dataArray->SetNumberOfComponents(1);
				dataArray->SetNumberOfTuples(opTriaHex.size());
				const Types::Real* cell_vec = opVec_Ele[i]->data();
				for (size_t j = 0; j < opTriaHex.size(); ++j) {
					dataArray->SetValue(j, cell_vec[j]);
				}
			}
			else if (opDataType_Ele[i] == DataOutType::vector) {
				// Vector data
				dataArray->SetNumberOfComponents(3); // Vector has 3 components
				dataArray->SetNumberOfTuples(opTriaHex.size());
				const Types::Real* cell_vec = opVec_Ele[i]->data();
				for (size_t j = 0; j < opTriaHex.size(); ++j) {
					dataArray->SetTuple3(j, cell_vec[3 * j], cell_vec[3 * j + 1], cell_vec[3 * j + 2]);
				}
			}
			else if (opDataType_Ele[i] == DataOutType::symmetric_tensor) {
				// Symmetric tensor data (3D has 6 independent components)
				dataArray->SetNumberOfComponents(6); // 6 components: σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃
				dataArray->SetNumberOfTuples(opTriaHex.size());
				const Types::Real* cell_vec = opVec_Ele[i]->data();
				for (size_t j = 0; j < opTriaHex.size(); ++j) {
					// Assuming data stored as [σ₁₁, σ₂₂, σ₃₃, σ₁₂, σ₁₃, σ₂₃]
					dataArray->SetTuple6(j, cell_vec[6 * j], cell_vec[6 * j + 1], cell_vec[6 * j + 2],
						cell_vec[6 * j + 3], cell_vec[6 * j + 4], cell_vec[6 * j + 5]);
				}
			}
			grid->GetCellData()->AddArray(dataArray);
		}

		// Add time field
		auto timeArray = vtkSmartPointer<vtkFloatArray>::New();
		timeArray->SetName("TimeValue");
		timeArray->SetNumberOfTuples(1);
		timeArray->SetValue(0, time);

		auto fieldData = vtkSmartPointer<vtkFieldData>::New();
		fieldData->AddArray(timeArray);
		grid->SetFieldData(fieldData);

		// Write to file
		auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName(filepath.c_str());  // Use full path
		writer->SetInputData(grid);
		writer->SetCompressorTypeToZLib();
		writer->Write();
	}

	void DataOut::write_vtu_virtualParticles(const unsigned int& result_number, const FSI_share_data& fsi_data) const
	{
		// 检查是否有数据
		if (fsi_data.FSI_virtualParticles_coordinates.empty()) {
			return; // 没有粒子，直接返回
		}

		// 确保输出目录存在
		struct stat info;
		const char* dirName = "ResultData";
		if (stat(dirName, &info) != 0 || !(info.st_mode & S_IFDIR)) {
#ifdef _WIN32
			if (_mkdir(dirName) != 0) {
				std::cerr << "Error: Could not create directory " << dirName << " (Windows)" << std::endl;
				return; // 创建失败则返回
			}
#else
			if (mkdir(dirName, 0755) != 0) {
				std::cerr << "Error: Could not create directory " << dirName << " (POSIX)" << std::endl;
				return; // 创建失败则返回
			}
#endif
		}

		// 文件名
		std::string fileName = std::string(dirName) + "/virParticles_" + std::to_string(result_number) + ".vtu";

		// 创建 VTK 数据结构
		auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		auto points = vtkSmartPointer<vtkPoints>::New();
		const auto& coordinates = fsi_data.FSI_virtualParticles_coordinates;
		vtkIdType nPoints = static_cast<vtkIdType>(coordinates.size());

		// --- 1. 定义点坐标 (使用 vtkDoubleArray 和 SetData) ---
		auto pointDataArray = vtkSmartPointer<vtkDoubleArray>::New();
		pointDataArray->SetNumberOfComponents(3); // X, Y, Z
		pointDataArray->Allocate(nPoints); // 预分配内存

		for (vtkIdType i = 0; i < nPoints; ++i) {
			// 检查 coordinates[i] 是否有效（可选，但更安全）
			// if (coordinates[i].size() == 3) { // std::array size is fixed at compile time
			pointDataArray->InsertNextTuple3(coordinates[i][0], coordinates[i][1], coordinates[i][2]);
			// } else {
			//	std::cerr << "Warning: Invalid coordinate data at index " << i << std::endl;
			//	pointDataArray->InsertNextTuple3(0.0, 0.0, 0.0); // Insert default value
			// }
		}
		// 将填充好的数据数组设置给 vtkPoints 对象
		points->SetData(pointDataArray);
		// 将 vtkPoints 对象设置给网格
		grid->SetPoints(points);
		// --- 结束点坐标设置 ---


		// --- 2. 定义单元 (顶点) ---
		auto cells = vtkSmartPointer<vtkCellArray>::New();
		cells->AllocateExact(nPoints, nPoints); // 预分配精确大小
		for (vtkIdType i = 0; i < nPoints; ++i) {
			vtkIdType pointId = i;
			cells->InsertNextCell(1, &pointId); // VTK_VERTEX 单元只有一个点
		}
		grid->SetCells(VTK_VERTEX, cells); // 设置单元类型和单元数组
		// --- 结束单元设置 ---


		// --- 3. 定义点数据 (速度、力等) ---
		const auto& velocities = fsi_data.FSI_virtualParticles_velocity;
		if (velocities.size() == nPoints) { // 检查大小匹配
			auto velocityArray = vtkSmartPointer<vtkFloatArray>::New(); // 可用 vtkDoubleArray
			velocityArray->SetName("velocity");
			velocityArray->SetNumberOfComponents(3);
			velocityArray->SetNumberOfTuples(nPoints);
			for (vtkIdType i = 0; i < nPoints; ++i) {
				velocityArray->SetTuple3(i, velocities[i][0], velocities[i][1], velocities[i][2]);
			}
			grid->GetPointData()->AddArray(velocityArray);
			grid->GetPointData()->SetActiveVectors("velocity"); // 标记为活动向量
		}
		else if (!velocities.empty()) { // 只有在非空但不匹配时才警告
			std::cerr << "Warning: Velocity data size (" << velocities.size()
				<< ") does not match coordinate size (" << nPoints
				<< ") in write_vtu_virtualParticles. Skipping velocity data." << std::endl;
		}

		const auto& forces = fsi_data.FSI_virtualParticles_nodeForce;
		if (forces.size() == nPoints) { // 检查大小匹配
			auto forceArray = vtkSmartPointer<vtkFloatArray>::New(); // 可用 vtkDoubleArray
			forceArray->SetName("force");
			forceArray->SetNumberOfComponents(3);
			forceArray->SetNumberOfTuples(nPoints);
			for (vtkIdType i = 0; i < nPoints; ++i) {
				forceArray->SetTuple3(i, forces[i][0], forces[i][1], forces[i][2]);
			}
			grid->GetPointData()->AddArray(forceArray);
		}
		else if (!forces.empty()) { // 只有在非空但不匹配时才警告
			std::cerr << "Warning: Force data size (" << forces.size()
				<< ") does not match coordinate size (" << nPoints
				<< ") in write_vtu_virtualParticles. Skipping force data." << std::endl;
		}
		// --- 结束点数据设置 ---


		// --- 4. 写入文件 ---
		auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName(fileName.c_str());
		writer->SetInputData(grid);
		writer->SetCompressorTypeToZLib(); // 启用压缩
		writer->SetDataModeToAppended();   // 使用附加数据模式

		if (!writer->Write()) {
			std::cerr << "Error writing VTU file: " << fileName << std::endl;
		}
		// --- 结束写入 ---
	}


}
