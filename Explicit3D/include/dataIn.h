#pragma once
#include <string>
#include <fstream>
#include <map>
#include "exDyna3D.h"


namespace EnSC {
	class exDyna3D;
	
	/**
	 * @brief 数据输入类，用于读取和解析输入文件
	 * 
	 * 该类负责从输入文件中读取模型数据并填充到exDyna3D类的数据结构中
	 */
	class DataIn {
	public:
		/**
		 * @brief 构造函数
		 * @param p_exdyna 引用到exDyna3D对象，用于存储读取的数据
		 */
		DataIn(exDyna3D& p_exdyna);
		
		/**
		 * @brief 读取输入文件
		 * @param fileName 输入文件名（.inp格式）
		 */
		void read_inp(std::string fileName);
	private:
		exDyna3D& exdyna; // 引用到求解器对象
		
		// 解析节点信息
		bool NODE();
		// 解析单元信息
		bool ELEMENT();
		// 解析材料属性
		bool MATERIAL();
		// 解析节点集合
		bool SET_NODE_LIST();
		// 解析边界条件
		bool BOUNDARY();
		// 解析初始条件
		bool INITIAL_CONDITIONS();
		// 解析分布载荷（包括重力）
		bool DLOAD();
		// 解析分布面载荷
		bool DSLOAD();
		// 解析时间步长信息
		bool TIME();
		// 解析幅值曲线
		bool AMP();
		// 解析单元集合
		bool SET_ELE_LIST();
		// 解析表面信息
		bool SURFACE();
		// 解析体积粘性参数
		bool BULK_VISCOSITY();
		// 解析节点位移约束
		void SPC_NODE();
		// 解析初始速度
		void INITIAL_VELOCITY();
		// 解析节点速度约束
		void VELOCITY_NODE();
		// 解析重力加载
		void GRAV(std::string& amp_name);
		
		std::string str; // 当前读取的行
		std::ifstream fin; // 输入文件流
		std::map<std::string, std::function<bool(DataIn*)>> k_func; // 关键字到处理函数的映射
	};
}