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
		 * @brief 解析状态枚举类，用于跟踪INP文件的解析状态
		 */
		enum class ParseState {
			GLOBAL,
			PART_DEFINITION,
			ASSEMBLY_DEFINITION,
			STEP_DEFINITION,
			MATERIAL_DEFINITION,
			PREDEFINED_FIELD,
			// 可以根据需要添加更多状态
		};

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
		ParseState currentState = ParseState::GLOBAL; // 当前解析状态
		std::string currentPartName;    // 当前Part名称
		std::string currentAssemblyName; // 当前Assembly名称
		std::string currentStepName;     // 当前Step名称
		std::string currentMaterialName; // 当前Material名称
		std::string currentInstanceName; // 当前实例名称
		std::string currentInstancePart; // 当前实例对应的部件
		
		// 节点和单元全局到局部索引的映射
		std::map<std::string, std::map<int, int>> instanceNodeMap; // 实例名 -> (全局索引 -> 局部索引)
		std::map<std::string, std::map<int, int>> instanceElemMap; // 实例名 -> (全局索引 -> 局部索引)
		
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
		// 解析输出间隔设置（自定义）
		bool OUTPUT_INTERVAL();
		// 解析Abaqus输出设置
		bool OUTPUT();
		// 解析节点位移约束
		void SPC_NODE();
		// 解析初始速度
		void INITIAL_VELOCITY();
		// 解析节点速度约束
		void VELOCITY_NODE();
		// 解析重力加载
		void GRAV(std::string& amp_name);
		// 解析部件实例
		bool INSTANCE();
		// 处理部件结束标记
		bool END_PART();
		// 处理装配结束标记
		bool END_ASSEMBLY();
		// 处理步骤结束标记
		bool END_STEP();
		// 处理实例结束标记
		bool END_INSTANCE();
		
		// 从关键字行提取名称（如NAME=value）
		std::string extractNameFromKeyword(const std::string& keywordLine, const std::string& prefix = "NAME=");
		// 解析包含实例引用的节点/单元集定义
		bool parseSetWithInstance(const std::string& str, std::string& setName, std::string& instanceName);
		
		std::string str; // 当前读取的行
		std::ifstream fin; // 输入文件流
		std::map<std::string, std::function<bool(DataIn*)>> k_func; // 关键字到处理函数的映射
	};
}