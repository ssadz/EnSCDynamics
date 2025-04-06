#include"../include/dataIn.h"
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include <string>
#include <cctype>

namespace EnSC {
	/**
	 * @brief 字符串转换为数值类型的模板函数声明
	 * 
	 * @tparam T 目标数值类型
	 * @param str 要转换的字符串
	 * @return 转换后的数值
	 */
	template<typename T>
	inline T convertString(const std::string& str);

	/**
	 * @brief 字符串转换为float的特化实现
	 * 
	 * @param str 要转换的字符串
	 * @return 转换后的float值
	 */
	template<>
	inline float convertString<float>(const std::string& str) {
		return std::stof(str);
	}

	/**
	 * @brief 字符串转换为double的特化实现
	 * 
	 * @param str 要转换的字符串
	 * @return 转换后的double值
	 */
	template<>
	inline double convertString<double>(const std::string& str) {
		return std::stod(str);
	}

	/**
	 * @brief 将字符串转换为大写
	 * 
	 * @param str 要转换的字符串（引用，将被修改）
	 */
	inline static void toUpperCase(std::string& str) {
		std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
	}

	/**
	 * @brief 移除字符串中的所有空白字符
	 * 
	 * @param str 要处理的字符串（引用，将被修改）
	 */
	inline static void removeSpaces(std::string& str) {
		str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
	}

	/**
	 * @brief 构造函数实现
	 * 
	 * 初始化关键字到处理函数的映射，用于解析输入文件中的不同部分
	 * 
	 * @param p_exdyna 引用到exDyna3D对象，用于存储解析结果
	 */
	DataIn::DataIn(exDyna3D& p_exdyna)
		:exdyna(p_exdyna) {
		// 建立关键字到处理函数的映射
		k_func["*NODE"] = &DataIn::NODE;                   // 节点数据
		k_func["*ELEMENT"] = &DataIn::ELEMENT;             // 单元数据
		k_func["*MATERIAL"] = &DataIn::MATERIAL;           // 材料属性
		k_func["*BOUNDARY"] = &DataIn::BOUNDARY;           // 边界条件
		k_func["*INITIAL CONDITION"] = &DataIn::INITIAL_CONDITIONS; // 初始条件
		k_func["*NSET"] = &DataIn::SET_NODE_LIST;          // 节点集合
		k_func["*DLOAD"] = &DataIn::DLOAD;                 // 分布载荷
		k_func["*DSLOAD"] = &DataIn::DSLOAD;               // 表面分布载荷
		k_func["*DYNAMIC"] = &DataIn::TIME;                // 动力学时间步长
		k_func["*AMPLITUDE"] = &DataIn::AMP;               // 幅值曲线
		k_func["*ELSET"] = &DataIn::SET_ELE_LIST;          // 单元集合
		k_func["*SURFACE"] = &DataIn::SURFACE;             // 表面定义
		k_func["*BULK VISCOSITY"] = &DataIn::BULK_VISCOSITY; // 体积粘性参数
	}

	/**
	 * @brief 读取输入文件实现
	 * 
	 * 打开并处理输入文件，根据关键字调用对应的处理函数
	 * 
	 * @param fileName 输入文件名
	 */
	void DataIn::read_inp(std::string fileName) {
		// 打开输入文件
		fin.open(fileName);
		if (!fin.is_open()) {
			std::cout << "无法打开文件！" << std::endl;
			exit(0);
		}
		else {
			std::cout << "成功打开文件 " << fileName << "！" << std::endl;
		}

		// 逐行读取并处理文件内容
		while (!fin.eof()) {
			if (fin.peek() == '*') {  // 检查是否是关键字行
				std::getline(fin, str);
				if (str.length() > 1) {
					char secondChar = str[1];
					if (secondChar != '*') {  // 忽略注释行（以 ** 开头）
						toUpperCase(str);  // 转换为大写，方便匹配
						for (const auto& pair : k_func) {
							const std::string& key = pair.first;
							// 检查 str 是否包含 key 作为子字符串
							if (str.find(key) != std::string::npos) {
								k_func[key](this);  // 调用对应的处理函数
								break;
							}
						}
					}
				}
			}
			else {
				std::getline(fin, str);  // 跳过非关键字行
			}
		}

		fin.close();
		std::cout << "成功读取文件 " << fileName << "！" << std::endl;
	}

	/**
	 * @brief 处理节点数据
	 * 
	 * 解析输入文件中的节点坐标数据
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::NODE() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*') {
			std::getline(fin, str);
		}

		// 逐行读取节点坐标
		Types::Point<3> x;
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) {
				continue;  // 跳过空行
			}

			std::istringstream iss(str);
			std::string token;

			// 读取节点ID（不使用）
			std::getline(iss, token, ',');

			// 读取X坐标
			std::getline(iss, token, ',');
			x[0] = convertString<Types::Real>(token);

			// 读取Y坐标
			std::getline(iss, token, ',');
			x[1] = convertString<Types::Real>(token);

			// 读取Z坐标
			std::getline(iss, token, ',');
			x[2] = convertString<Types::Real>(token);

			// 将节点坐标添加到节点数组
			exdyna.vertices.emplace_back(x);
		}
		return true;
	}

	/**
	 * @brief 处理材料属性数据
	 * 
	 * 解析输入文件中的材料属性数据，包括密度、弹性模量和泊松比
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::MATERIAL() {
		Types::Real density = 0.0;          // 材料密度
		Types::Real elasticModulus = 0.0;   // 弹性模量
		Types::Real poissonRatio = 0.0;     // 泊松比
		std::string materialName;           // 材料名称

		// 局部定义的大写转换函数
		auto toUpperCase = [](std::string& str) {
			std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
			};

		// 读取并解析材料块
		while (std::getline(fin, str)) {
			toUpperCase(str);
			if (str.find("*MATERIAL") != std::string::npos) {
				// 获取材料名称
				auto pos = str.find("NAME=");
				if (pos != std::string::npos) {
					materialName = str.substr(pos + 5);
				}
			}
			else if (str.find("*DENSITY") != std::string::npos) {
				// 读取密度
				std::getline(fin, str);
				std::istringstream densityStream(str);
				densityStream >> density;
			}
			else if (str.find("*ELASTIC") != std::string::npos) {
				// 读取弹性参数（弹性模量和泊松比）
				std::getline(fin, str);
				std::istringstream iss(str);
				std::string token;
				std::getline(iss, token, ',');
				elasticModulus = convertString<Types::Real>(token);
				std::getline(iss, token, ',');
				poissonRatio = convertString<Types::Real>(token);
			}
			else if (str.find("*") == 0 && str.find("*MATERIAL") == std::string::npos) {
				// 遇到其他关键字，表示材料块结束
				break;
			}
		}

		// 将材料属性设置到模型中
		exdyna.mMatElastic.rho = density;
		exdyna.mMatElastic.E = elasticModulus;
		exdyna.mMatElastic.v = poissonRatio;
		exdyna.mMatElastic.update();  // 更新派生属性（如剪切模量、体积模量等）

		return true;
	}

	/**
	 * @brief 处理单元数据
	 * 
	 * 解析输入文件中的单元定义数据，包括单元类型和节点索引
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::ELEMENT() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*') {
			std::getline(fin, str);
		}

		int PID = 1;  // 默认物理ID
		int MID = 1;  // 默认材料ID
		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);

			// 跳过空行
			if (str.empty()) continue;

			std::istringstream iss(str);
			std::string token;

			// 读取单元ID（不使用）
			std::getline(iss, token, ',');

			// 读取单元的节点索引
			std::vector<Types::Vertex_index> verticesIndex;
			while (std::getline(iss, token, ',')) {
				verticesIndex.push_back(std::stoul(token));
			}

			// 节点索引从1开始，转换为从0开始
			for (auto& idx : verticesIndex) {
				idx -= 1;
			}

			// 创建新单元并设置属性
			exdyna.hexahedron_elements.emplace_back();
			exdyna.hexahedron_elements.back().set_PID(PID);
			exdyna.hexahedron_elements.back().set_MID(MID);
			exdyna.hexahedron_elements.back().set_verticesIndex(verticesIndex);
		}
		return true;
	}

	/**
	 * @brief 处理节点集合数据
	 * 
	 * 解析输入文件中的节点集合定义，用于后续引用特定的节点组
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::SET_NODE_LIST() {
		// 从关键字行中提取节点集合名称
		std::string::size_type pos = str.find("NSET=");
		pos += 5;
		std::string set_name;
		while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
			set_name += str[pos++];
		}

		bool is_generate = false;
		if (str.find("GENERATE") != std::string::npos) {
			is_generate = true;
		}

		// 读取节点集合数据
		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);
			std::vector<std::size_t> set_nodes;
			std::istringstream iss(str);
			while (!iss.eof()) {
				std::size_t value;
				iss >> value;
				set_nodes.push_back(value);
				if (iss.peek() == ',') {
					iss.ignore();
				}
			}

			// 处理生成模式（start,end,step）
			if (is_generate && set_nodes.size() == 3) {
				std::size_t start = set_nodes[0];
				std::size_t end = set_nodes[1];
				std::size_t step = set_nodes[2];
				set_nodes.clear();
				for (std::size_t i = start; i <= end; i += step) {
					set_nodes.push_back(i);
				}
			}

			// 索引从1开始，转换为从0开始
			for (auto& idx : set_nodes) {
				idx -= 1;
			}

			// 存储节点集合
			exdyna.map_set_node_list[set_name] = set_nodes;
		}
		return true;
	}

	bool DataIn::BOUNDARY() {
		if (str.find("TYPE") != std::string::npos) {
			if (str.find("VELOCITY") != std::string::npos) {
				VELOCITY_NODE();
			}
		}
		else {
			SPC_NODE();
		}
		return true;
	}

	/**
	 * @brief 处理初始条件数据
	 * 
	 * 解析输入文件中的初始条件数据，根据类型调用相应的处理函数
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::INITIAL_CONDITIONS() {
		// 检查初始条件类型并调用相应的处理函数
		if (str.find("TYPE") != std::string::npos) {
			if (str.find("VELOCITY") != std::string::npos) {
				// 处理初始速度条件
				INITIAL_VELOCITY();
			}
		}
		return true;
	}

	/**
	 * @brief 处理分布载荷数据
	 * 
	 * 解析输入文件中的分布载荷数据，包括重力等
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::DLOAD() {
		// 从关键字行中提取幅值曲线名称
		std::string amp_name;
		std::string::size_type pos = str.find("AMPLITUDE=");
		if (pos != std::string::npos) {
			pos += 10;
			while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
				amp_name += str[pos++];
			}
		}
		
		// 读取分布载荷数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.find("GRAV") != std::string::npos) {
				// 处理重力载荷
				GRAV(amp_name);
			}
		}

		return true;
	}

	/**
	 * @brief 处理动力学时间步长数据
	 * 
	 * 解析输入文件中的动力学时间步长参数，设置总模拟时间
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::TIME() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*') {
			std::getline(fin, str);
		}

		// 读取时间步长数据
		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);
			std::istringstream iss(str);
			std::string token;
			
			// 跳过第一个值
			std::getline(iss, token, ',');
			
			// 读取总模拟时间
			std::getline(iss, token, ',');
			exdyna.totalTime = convertString<Types::Real>(token);
		}
		return true;
	}

	/**
	 * @brief 处理体积粘性参数数据
	 * 
	 * 解析输入文件中的体积粘性参数，用于控制冲击波传播
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::BULK_VISCOSITY() {
		std::getline(fin, str);
		std::istringstream iss(str);
		std::string token;

		// 读取线性体积粘性系数
		std::getline(iss, token, ',');
		exdyna.Cvisl = convertString<Types::Real>(token);
		
		// 读取二次体积粘性系数
		std::getline(iss, token, ',');
		exdyna.Cvisq = convertString<Types::Real>(token);
		
		return true;
	}

	/**
	 * @brief 处理单点约束节点数据
	 * 
	 * 解析输入文件中的单点约束节点数据，用于固定节点位移
	 */
	void DataIn::SPC_NODE() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*')
			std::getline(fin, str);

		// 读取单点约束数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			exdyna.boundary_spc_node.emplace_back();
			std::istringstream iss(str);
			std::string token;

			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				// 使用预定义的节点集合
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().first = exdyna.map_set_node_list[token];
			}
			else {
				// 单个节点
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().first.push_back(std::stoi(token));
			}

			// 处理特殊约束类型（固定所有自由度）
			if (str.find("PINNED") != std::string::npos || str.find("ENCASTRE") != std::string::npos) {
				exdyna.boundary_spc_node.back().second.at(0) = 0;  // 起始自由度
				exdyna.boundary_spc_node.back().second.at(1) = 2;  // 结束自由度
				exdyna.boundary_spc_node.back().second.at(2) = 1;  // 增量
			}
			else {
				// 读取自由度约束定义
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().second.at(0) = std::stoi(token) - 1;  // 起始自由度
				
				std::getline(iss, token, ',');
				if (token.empty()) {
					// 如果没有指定结束自由度，则设置为与起始自由度相同
					exdyna.boundary_spc_node.back().second.at(1) = exdyna.boundary_spc_node.back().second.at(0);
				}
				else {
					exdyna.boundary_spc_node.back().second.at(1) = std::stoi(token) - 1;  // 结束自由度
				}
				
				// 读取增量
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().second.at(2) = std::stoi(token);
			}
		}
	}

	/**
	 * @brief 处理初始速度数据
	 * 
	 * 解析输入文件中的初始速度数据，为节点设置初始速度
	 */
	void DataIn::INITIAL_VELOCITY() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*')
			std::getline(fin, str);
		
		// 读取初始速度数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			exdyna.ini_vel_generation.emplace_back();
			std::istringstream iss(str);
			std::string token;

			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				// 使用预定义的节点集合
				std::getline(iss, token, ',');
				exdyna.ini_vel_generation.back().first = exdyna.map_set_node_list[token];
			}
			else {
				// 单个节点
				std::getline(iss, token, ',');
				exdyna.ini_vel_generation.back().first.push_back(std::stoi(token) - 1);
			}

			// 读取自由度索引和速度值
			std::getline(iss, token, ',');
			exdyna.ini_vel_generation.back().second.first = std::stoi(token) - 1;  // 自由度索引
			std::getline(iss, token, ',');
			exdyna.ini_vel_generation.back().second.second = convertString<Types::Real>(token);  // 速度值
		}
	}

	/**
	 * @brief 处理速度边界条件节点数据
	 * 
	 * 解析输入文件中的速度边界条件数据，为节点设置速度约束
	 */
	void DataIn::VELOCITY_NODE() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*')
			std::getline(fin, str);
		
		// 读取速度边界条件数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			std::istringstream iss(str);
			std::string token;

			// 创建临时结构存储速度边界条件
			std::pair<std::array<std::size_t, 2>, Types::Real>temp_pair;
			temp_pair.first[0] = 0;
			temp_pair.first[1] = 0;
			temp_pair.second = 0.0;
			
			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				std::getline(iss, token, ',');
				exdyna.boundary_vel_node.emplace_back(exdyna.map_set_node_list[token], temp_pair);
			}
			else {
				std::getline(iss, token, ',');
				std::vector<std::size_t>temp_vector;
				temp_vector.emplace_back(std::stoi(token) - 1);
				exdyna.boundary_vel_node.emplace_back(temp_vector, temp_pair);
			}

			std::getline(iss, token, ',');
			exdyna.boundary_vel_node.back().second.first[0] = std::stoi(token) - 1;
			std::getline(iss, token, ',');
			if (token.empty()) {
				exdyna.boundary_vel_node.back().second.first[1] = exdyna.boundary_vel_node.back().second.first[0];
			}
			else {
				exdyna.boundary_vel_node.back().second.first[1] = std::stoi(token) - 1;
			}

			if (std::getline(iss, token, ',')) {
				exdyna.boundary_vel_node.back().second.second = convertString<Types::Real>(token);
			}
		}
	}

	//要修改
	void DataIn::GRAV(std::string& amp_name) {
		std::istringstream iss(str);
		std::string token;
		std::getline(iss, token, ',');
		std::getline(iss, token, ',');
		std::getline(iss, token, ',');

		Types::Real value = convertString<Types::Real>(token);
		std::getline(iss, token, ',');
		Types::Real direction_x = convertString<Types::Real>(token);
		std::getline(iss, token, ',');
		Types::Real direction_y = convertString<Types::Real>(token);
		std::getline(iss, token, ',');
		Types::Real direction_z = convertString<Types::Real>(token);

		std::get<0>(exdyna.gravity) = true;
		std::get<1>(exdyna.gravity) = amp_name;
		std::get<2>(exdyna.gravity) = value;
		std::get<3>(exdyna.gravity) = direction_x;
		std::get<4>(exdyna.gravity) = direction_y;
		std::get<5>(exdyna.gravity) = direction_z;

	}

	bool DataIn::AMP() {
		std::string::size_type pos = str.find("NAME=");
		pos += 5;
		std::string amp_name;
		while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
			amp_name += str[pos++];
		}

		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);
			std::istringstream iss(str);
			std::istringstream iss_c(str);
			std::string token;

			int N_Real = 0;
			while (std::getline(iss, token, ',')) {
				N_Real++;
			}

			std::vector<std::array<Types::Real, 2>>temp;
			N_Real /= 2;
			temp.resize(N_Real);

			for (int i = 0; i < temp.size(); i++) {
				std::getline(iss_c, token, ',');
				temp[i][0] = convertString<Types::Real>(token);
				std::getline(iss_c, token, ',');
				temp[i][1] = convertString<Types::Real>(token);
			}
			exdyna.map_amp_list[amp_name] = temp;
		}

		return true;
	}

	bool DataIn::SET_ELE_LIST() {
		std::string::size_type pos = str.find("ELSET=");
		pos += 6;
		std::string set_ele_name;
		while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
			set_ele_name += str[pos++];
		}

		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);


			std::istringstream iss(str);
			std::istringstream iss_c(str);
			std::string token;
			std::vector<std::size_t>temp_index_set;
			if (str.find("EALL") != std::string::npos) {
				temp_index_set.clear();
				for (std::size_t i = 1; i <= exdyna.hexahedron_elements.size(); i++) {
					temp_index_set.push_back(i);
				}
			}
			else {
				int N_Real = 0;
				while (std::getline(iss, token, ',')) {
					N_Real++;
				}

				if (N_Real == 1 && 2) {
					std::getline(iss_c, token, ',');
					std::size_t index = std::stoi(token);

					temp_index_set.push_back(index);
					if (N_Real == 2) {
						std::getline(iss_c, token, ',');
						index = std::stoi(token);
						temp_index_set.push_back(index);
					}
				}

				if (N_Real == 3) {
					std::getline(iss_c, token, ',');
					std::size_t index = std::stoi(token);
					temp_index_set.push_back(index);
					std::getline(iss_c, token, ',');
					index = std::stoi(token);
					temp_index_set.push_back(index);
					std::getline(iss_c, token, ',');
					index = std::stoi(token);
					temp_index_set.push_back(index);

					std::size_t start = temp_index_set[0];
					std::size_t end = temp_index_set[1];
					std::size_t step = temp_index_set[2];
					temp_index_set.clear();
					for (std::size_t i = start; i <= end; i += step) {
						temp_index_set.push_back(i);
					}
				}
			}
			for (auto& idx : temp_index_set) {
				idx -= 1;
			}
			exdyna.map_set_ele_list[set_ele_name] = temp_index_set;
		}
		return true;
	}

	bool DataIn::SURFACE() {
		std::string::size_type pos = str.find("NAME=");
		if (pos != std::string::npos) {
			pos += 5;
			std::string surf_name;
			while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
				surf_name += str[pos++];
			}
			// 移除 surf_name 中的空格
			removeSpaces(surf_name);

			while (fin.peek() != '*' && fin.peek() != EOF) {
				std::getline(fin, str);
				toUpperCase(str);
				std::istringstream iss(str);
				std::string token;

				std::getline(iss, token, ',');
				std::string ele_set_name = token;
				// 移除 ele_set_name 中的空格
				removeSpaces(ele_set_name);

				std::getline(iss, token, ',');
				std::string surf_type = token;
				// 移除 surf_type 中的空格
				removeSpaces(surf_type);

				std::pair<std::string, std::string> temp_surf_name_type(ele_set_name, surf_type);
				exdyna.map_set_surface_list[surf_name] = temp_surf_name_type;
			}
		}
		return true;
	}


	bool DataIn::DSLOAD() {
		std::string amp_name;
		std::string::size_type pos = str.find("AMPLITUDE=");
		if (pos != std::string::npos) {
			pos += 10;
			while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
				amp_name += str[pos++];
			}
		}

		// Remove spaces from amp_name
		removeSpaces(amp_name);

		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);
			toUpperCase(str);
			std::istringstream iss(str);
			std::string token;

			std::getline(iss, token, ',');
			std::string surf_name = token;
			removeSpaces(surf_name);

			std::getline(iss, token, ',');
			std::string load_type = token;
			removeSpaces(load_type);

			std::getline(iss, token, ',');
			Types::Real load_value = convertString<Types::Real>(token);
			exdyna.dsload.emplace_back();

			std::get<0>(exdyna.dsload.back()) = surf_name;
			std::get<1>(exdyna.dsload.back()) = load_type;
			std::get<2>(exdyna.dsload.back()) = amp_name;
			std::get<3>(exdyna.dsload.back()) = load_value;
		}
		return true;
	}

}

