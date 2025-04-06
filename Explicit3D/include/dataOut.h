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
		/**
		 * @brief 定义输出数据类型枚举
		 */
		enum DataType { 
			scalar,           // 标量数据
			vector,           // 向量数据
			symmetric_tensor  // 对称张量数据
		};
	}

	/**
	 * @brief 数据输出类，负责将计算结果输出为可视化文件
	 */
	class DataOut {
	public:
		/**
		 * @brief 构造函数
		 * @param pVerticesAll 节点坐标数组
		 * @param pTriaHex 单元数组
		 * @param time 当前时间
		 */
		DataOut(const Types::VerticesAll<3>& pVerticesAll, const std::vector<Element_HexN8>& pTriaHex, Types::Real time);
		~DataOut() {}

		/**
		 * @brief 添加节点数据
		 * @param dataType 数据类型（标量、向量或对称张量）
		 * @param vec 数据数组指针
		 * @param name 数据名称
		 */
		void add_node_data(DataOutType::DataType dataType, const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>* vec, std::string name);

		/**
		 * @brief 添加单元数据
		 * @param dataType 数据类型（标量、向量或对称张量）
		 * @param vec 数据数组指针
		 * @param name 数据名称
		 */
		void add_ele_data(DataOutType::DataType dataType, const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>* vec, std::string name);
		
		/**
		 * @brief 将数据写入VTU文件（VTK非结构网格格式）
		 * @param fileName 输出文件名
		 */
		void write_vtu(const std::string& fileName);
		
		/**
		 * @brief 将虚拟粒子数据写入VTU文件，用于FSI可视化
		 * @param result_number 结果编号
		 * @param fsi_data FSI共享数据
		 */
		void write_vtu_virtualParticles(const unsigned int& result_number, const FSI_share_data& fsi_data) const;
	private:
		const Types::VerticesAll<3>& opVerticesAll;    // 节点坐标数组引用
		const std::vector<Element_HexN8>& opTriaHex;   // 单元数组引用
		Types::Real time;                              // 当前时间

		// 节点数据存储
		std::vector<DataOutType::DataType> opDataType;  // 节点数据类型
		std::vector<const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>*> opVec;  // 节点数据指针
		std::vector<std::string> opName;                // 节点数据名称

		// 单元数据存储
		std::vector<DataOutType::DataType> opDataType_Ele;  // 单元数据类型
		std::vector<const Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>*> opVec_Ele;  // 单元数据指针
		std::vector<std::string> opName_Ele;               // 单元数据名称
	};

}
#endif // DATAOUT_H
