#include"../include/dataIn.h"
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include <string>
#include <cctype>
#include "../include/InpData.h"

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
		:exdyna(p_exdyna), currentPartPtr(nullptr), currentAssemblyPtr(nullptr) {
		// 初始化关键字映射表，将输入文件中的关键字与对应的处理函数关联
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
		k_func["*OUTPUT INTERVAL"] = &DataIn::OUTPUT_INTERVAL; // 自定义输出间隔设置
		k_func["*OUTPUT"] = &DataIn::OUTPUT;               // Abaqus输出设置
		k_func["*INSTANCE"] = &DataIn::INSTANCE;           // 部件实例
		k_func["*END PART"] = &DataIn::END_PART;           // 部件结束标记
		k_func["*END ASSEMBLY"] = &DataIn::END_ASSEMBLY;   // 装配结束标记
		k_func["*END STEP"] = &DataIn::END_STEP;           // 步骤结束标记
		k_func["*END INSTANCE"] = &DataIn::END_INSTANCE;   // 实例结束标记
	}

	/**
	 * @brief 从关键字行提取名称
	 * 
	 * @param keywordLine 关键字行
	 * @param prefix 名称前缀（如"NAME="）
	 * @return 提取的名称，如果未找到则返回空字符串
	 */
	std::string DataIn::extractNameFromKeyword(const std::string& keywordLine, const std::string& prefix) {
		std::string result = "";
		size_t namePos = keywordLine.find(prefix);
		if (namePos != std::string::npos) {
			namePos += prefix.length();
			size_t endPos = keywordLine.find(",", namePos);
			if (endPos == std::string::npos) {
				endPos = keywordLine.length();
			}
			result = keywordLine.substr(namePos, endPos - namePos);
			// 去除可能的空格
			result.erase(std::remove_if(result.begin(), result.end(), ::isspace), result.end());
		}
		return result;
	}

	/**
	 * @brief 解析包含实例引用的节点/单元集定义
	 * 
	 * @param str 关键字行，包含instance=xxx的定义
	 * @param setName 输出参数：集合名称
	 * @param instanceName 输出参数：实例名称
	 * @return bool 解析是否成功
	 */
	bool DataIn::parseSetWithInstance(const std::string& str, std::string& setName, std::string& instanceName) {
		// 提取集合名称
		setName = extractNameFromKeyword(str, "NSET=");
		if (setName.empty()) {
			setName = extractNameFromKeyword(str, "ELSET=");
		}
		
		// 提取实例名称
		instanceName = extractNameFromKeyword(str, "INSTANCE=");
		
		return !setName.empty();
	}

	/**
	 * @brief 处理部件实例
	 * 
	 * 解析装配中的部件实例定义
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::INSTANCE() {
		if (!currentAssemblyPtr) return false;
		std::string instance_name = extractNameFromKeyword(str, "NAME=");
		std::string part_name = extractNameFromKeyword(str, "PART=");
		InpInstance inst;
		inst.name = instance_name;
		inst.part_name = part_name;
		
		// 初始化新增的成员
		inst.translation = {0.0, 0.0, 0.0}; // 默认无平移
		inst.node_start_index = 0;          // 将在transferToExDyna中设置
		inst.element_start_index = 0;       // 将在transferToExDyna中设置
		
		// 处理可能的平移信息（这里简化处理，实际中可能需要解析更多参数）
		// TODO: 实现解析ABAQUS INP文件中的平移和旋转变换
		
		currentAssemblyPtr->instances.push_back(inst);
		return true;
	}

	/**
	 * @brief 处理部件结束标记
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::END_PART() {
		currentPartPtr = nullptr;
		currentState = ParseState::GLOBAL;
		std::cout << "结束部件定义: " << currentPartName << std::endl;
		return true;
	}

	/**
	 * @brief 处理装配结束标记
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::END_ASSEMBLY() {
		currentAssemblyPtr = nullptr;
		currentState = ParseState::GLOBAL;
		std::cout << "结束装配定义: " << currentAssemblyName << std::endl;
		return true;
	}

	/**
	 * @brief 处理步骤结束标记
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::END_STEP() {
		currentState = ParseState::GLOBAL;
		std::cout << "结束步骤定义: " << currentStepName << std::endl;
		
		// 增加更多调试信息
		if (!exdyna.steps.empty()) {
			StepData& currentStep = exdyna.steps.back();
			std::cout << "------------------步骤调试信息------------------" << std::endl;
			std::cout << "步骤名称: " << currentStep.name << std::endl;
			std::cout << "时间周期: " << currentStep.timePeriod << std::endl;
			std::cout << "本步骤定义的固定约束节点条件数: " << currentStep.boundary.spc_nodes.size() << std::endl;
			std::cout << "本步骤定义的速度约束节点条件数: " << currentStep.boundary.vel_nodes.size() << std::endl;
			
			// 打印重置标记信息
			if (currentStep.resetSpcBoundary) {
				std::cout << "本步骤重置了所有位移边界条件 (Boundary, op=NEW)" << std::endl;
			}
			if (currentStep.resetVelBoundary) {
				std::cout << "本步骤重置了所有速度边界条件 (Boundary, op=NEW, type=VELOCITY)" << std::endl;
			}
			if (currentStep.resetDload) {
				std::cout << "本步骤重置了所有分布载荷 (DLOAD, op=NEW)" << std::endl;
			}
			if (currentStep.resetDsload) {
				std::cout << "本步骤重置了所有分布面载荷 (DSLOAD, op=NEW)" << std::endl;
			}
			
			// 提示关于边界条件继承的情况
			if (exdyna.steps.size() > 1) {
				if (!currentStep.resetSpcBoundary && !currentStep.resetVelBoundary &&
					!currentStep.resetDload && !currentStep.resetDsload) {
					std::cout << "注意: 在运行时，如果本步骤未定义某类条件，将从前一步骤继承" << std::endl;
				} else {
					std::cout << "注意: 本步骤有重置标记，部分条件将不会从前面步骤继承" << std::endl;
				}
			}
			
			std::cout << "------------------------------------------------" << std::endl;
		} else {
			std::cout << "警告: 当前步骤数据未保存!" << std::endl;
		}
		
		return true;
	}

	/**
	 * @brief 处理实例结束标记
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::END_INSTANCE() {
		currentInstanceName = "";
		currentInstancePart = "";
		std::cout << "结束实例定义" << std::endl;
		return true;
	}

	/**
	 * @brief 读取输入文件实现
	 * 
	 * 打开并处理输入文件，根据关键字调用对应的处理函数
	 * 优化为更健壮的解析方式，与Abaqus INP文件结构对应
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

		// 初始化解析状态和名称
		currentState = ParseState::GLOBAL;
		currentPartName = "";
		currentAssemblyName = "";
		currentStepName = "";
		currentMaterialName = "";
		currentInstanceName = "";
		currentInstancePart = "";
		currentPartPtr = nullptr;
		currentAssemblyPtr = nullptr;

		// 逐行读取并处理文件内容
		while (!fin.eof()) {
			if (fin.peek() == '*') {  // 检查是否是关键字行
				std::getline(fin, str);
				if (str.length() > 1) {
					char secondChar = str[1];
					if (secondChar != '*') {  // 忽略注释行（以 ** 开头）
						toUpperCase(str);  // 转换为大写，方便匹配
						std::cout << "处理关键字行: " << str << std::endl;  // 添加调试输出
						
						// 用于控制循环流程的变量
						bool continue_main_loop = false;
						
						// 判断主要结构关键字并更新状态
						if (str.find("*HEADING") != std::string::npos) {
							// 文件头部，可以跳过或处理特定信息
							while (fin.peek() != '*' && !fin.eof()) {
								std::getline(fin, str);  // 跳过heading内容
							}
							continue;
						}
						else if (str.find("*PART") != std::string::npos) {
							currentPartName = extractNameFromKeyword(str);
							inp_data.parts.emplace_back();
							currentPartPtr = &inp_data.parts.back();
							currentPartPtr->name = currentPartName;
							currentState = ParseState::PART_DEFINITION;
							std::cout << "进入部件定义模式: " << currentPartName << std::endl;
							continue;
						}
						else if (str.find("*ASSEMBLY") != std::string::npos) {
							currentAssemblyName = extractNameFromKeyword(str);
							if (currentAssemblyName.empty()) {
								currentAssemblyName = "MainAssembly"; // 默认名称
							}
							inp_data.assemblies.emplace_back();
							currentAssemblyPtr = &inp_data.assemblies.back();
							currentAssemblyPtr->name = currentAssemblyName;
							currentState = ParseState::ASSEMBLY_DEFINITION;
							std::cout << "进入装配定义模式: " << currentAssemblyName << std::endl;
							continue;
						}
						else if (str.find("*MATERIAL") != std::string::npos) {
							currentMaterialName = extractNameFromKeyword(str);
							currentState = ParseState::MATERIAL_DEFINITION;
							std::cout << "进入材料定义模式: " << currentMaterialName << std::endl;
							continue;
						}
						else if (str.find("*STEP") != std::string::npos) {
							currentStepName = extractNameFromKeyword(str);
							if (currentStepName.empty()) {
								currentStepName = "Step-" + std::to_string(exdyna.steps.size() + 1); // 默认名称
							}
							currentState = ParseState::STEP_DEFINITION;
							// 创建新的步骤数据
							StepData stepData;
							stepData.name = currentStepName;
							
							// --- 不需要在这里显式继承边界条件，切换步骤时会自动处理 ---
							// 每个步骤只保存它自己显式定义的边界条件
							// 继承逻辑在 setCurrentStep 方法中处理
							
							exdyna.steps.push_back(stepData);
							std::cout << "进入步骤定义模式: " << currentStepName << std::endl;
							continue;
						}
						else if (str.find("*INITIAL CONDITIONS") != std::string::npos) {
							currentState = ParseState::PREDEFINED_FIELD;
							std::cout << "进入预定义场模式" << std::endl;
						}
						
						// 直接检查是否有对应的处理函数，包括结束标记和实例定义
						bool handled = false;
						for (const auto& pair : k_func) {
							const std::string& key = pair.first;
							if (str.find(key) != std::string::npos) {
								k_func[key](this);
								handled = true;
								break;
							}
						}
						
						// 如果没有找到直接对应的处理函数，则根据当前状态分类处理
						if (!handled) {
							switch (currentState) {
							case ParseState::GLOBAL:
								// 全局上下文，无特殊处理
								break;
								
							case ParseState::PART_DEFINITION:
								// 在部件定义上下文中处理特定关键字
								if (str.find("*SOLID SECTION") != std::string::npos) {
									// 处理实体截面定义
									std::string elset = extractNameFromKeyword(str, "ELSET=");
									std::string material = extractNameFromKeyword(str, "MATERIAL=");
									std::cout << "为单元集 " << elset << " 应用材料 " << material << std::endl;
									
									// 跳过截面定义行
									if (fin.peek() != '*' && !fin.eof()) {
										std::getline(fin, str);
									}
								}
								break;
								
							case ParseState::MATERIAL_DEFINITION:
								// 在材料定义上下文中处理特定关键字
								if (str.find("*DENSITY") != std::string::npos) {
									// 读取密度值
									if (fin.peek() != '*' && !fin.eof()) {
										std::getline(fin, str);
										std::istringstream iss(str);
										Types::Real density;
										iss >> density;
										std::cout << "设置材料 " << currentMaterialName << " 的密度: " << density << std::endl;
										exdyna.mMatElastic.rho = density;
									}
								}
								else if (str.find("*ELASTIC") != std::string::npos) {
									// 读取弹性参数
									if (fin.peek() != '*' && !fin.eof()) {
										std::getline(fin, str);
										std::istringstream iss(str);
										std::string token;
										std::getline(iss, token, ',');
										Types::Real elasticModulus = convertString<Types::Real>(token);
										std::getline(iss, token, ',');
										Types::Real poissonRatio = convertString<Types::Real>(token);
										
										std::cout << "设置材料 " << currentMaterialName 
												  << " 的弹性模量: " << elasticModulus 
												  << "，泊松比: " << poissonRatio << std::endl;
										
										exdyna.mMatElastic.E = elasticModulus;
										exdyna.mMatElastic.v = poissonRatio;
										exdyna.mMatElastic.update();  // 更新派生属性
									}
								}
								break;
								
							case ParseState::ASSEMBLY_DEFINITION:
								// 在装配定义上下文中，处理实例中的节点集和单元集
								// 这里可能需要处理instance引用
								break;
								
							case ParseState::STEP_DEFINITION:
								// 特殊处理*BOUNDARY命令，以确保连续的边界命令能够被单独处理
								if (str.find("*BOUNDARY") != std::string::npos) {
									std::cout << "**特殊处理步骤内边界条件命令: " << str << std::endl;
									BOUNDARY();  // 直接调用处理函数
									std::cout << "**步骤内边界条件命令处理完成" << std::endl;
									// 继续外部循环，跳过其他处理
									handled = true;  // 标记为已处理
									continue_main_loop = true;  // 设置一个标志，指示应该继续主循环
									break;  // 跳出switch
								}
								
								// 在步骤定义上下文中，处理输出请求等
								if (str.find("*RESTART") != std::string::npos) {
									// 处理重启输出设置
									while (fin.peek() != '*' && !fin.eof()) {
										std::getline(fin, str);  // 跳过restart内容
									}
								}
								else if (str.find("*OUTPUT, FIELD") != std::string::npos ||
										 str.find("*OUTPUT, HISTORY") != std::string::npos) {
									// 处理场输出或历史输出设置
									if (str.find("NUMBER INTERVAL=") != std::string::npos) {
										std::string intervalStr = extractNameFromKeyword(str, "NUMBER INTERVAL=");
										if (!intervalStr.empty()) {
											int interval = std::stoi(intervalStr);
											exdyna.time_interval = exdyna.totalTime / interval;
											exdyna.time_interval_set = true;
											std::cout << "从输出设置获取间隔数: " << interval 
													  << "，设置time_interval=" << exdyna.time_interval << std::endl;
										}
									}
									
									// 跳过输出变量定义
									while (fin.peek() != '*' && !fin.eof()) {
										std::getline(fin, str);
									}
								}
								break;
								
							case ParseState::PREDEFINED_FIELD:
								// 处理预定义场，如初始条件
								if (str.find("*INITIAL CONDITIONS") != std::string::npos) {
									if (str.find("TYPE=VELOCITY") != std::string::npos) {
										INITIAL_VELOCITY();
									}
								}
								break;
							}
						}
						
						// 如果设置了continue_main_loop标志，则跳过后续处理，直接进入下一次循环
						if (continue_main_loop) {
							continue;
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
		
		// 新增：调用transferToExDyna()将inp_data转换为exdyna数据
		transferToExDyna();
		
		// 更新材料属性，确保在所有数据读取后执行
		exdyna.mMatElastic.update();
		
		// 打印解析统计信息
		std::cout << "解析完成: " << std::endl;
		if (!currentPartName.empty()) {
			std::cout << "- 部件名称: " << currentPartName << std::endl;
		}
		if (!currentAssemblyName.empty()) {
			std::cout << "- 装配名称: " << currentAssemblyName << std::endl;
		}
		if (!currentStepName.empty()) {
			std::cout << "- 步骤名称: " << currentStepName << std::endl;
		}
		if (!currentMaterialName.empty()) {
			std::cout << "- 材料名称: " << currentMaterialName << std::endl;
		}
		std::cout << "- 节点数量: " << exdyna.vertices.size() << std::endl;
		std::cout << "- 单元数量: " << exdyna.hexahedron_elements.size() << std::endl;
	}

	/**
	 * @brief 处理节点数据
	 * 
	 * 解析输入文件中的节点坐标数据
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::NODE() {
		if (!currentPartPtr) return false;
		while (fin.peek() == '*') std::getline(fin, str);
		Types::Point<3> x;
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) continue;
			std::istringstream iss(str);
			std::string token;
			std::getline(iss, token, ','); // 节点ID
			std::getline(iss, token, ','); x[0] = convertString<Types::Real>(token);
			std::getline(iss, token, ','); x[1] = convertString<Types::Real>(token);
			std::getline(iss, token, ','); x[2] = convertString<Types::Real>(token);
			currentPartPtr->nodes.push_back(x);
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
		InpMaterial material;
		material.name = currentMaterialName;
		
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
					material.name = str.substr(pos + 5);
				}
			}
			else if (str.find("*DENSITY") != std::string::npos) {
				// 读取密度
				std::getline(fin, str);
				std::istringstream densityStream(str);
				densityStream >> material.density;
			}
			else if (str.find("*ELASTIC") != std::string::npos) {
				// 读取弹性参数（弹性模量和泊松比）
				std::getline(fin, str);
				std::istringstream iss(str);
				std::string token;
				std::getline(iss, token, ',');
				material.E = convertString<Types::Real>(token);
				std::getline(iss, token, ',');
				material.v = convertString<Types::Real>(token);
			}
			else if (str.find("*") == 0 && str.find("*MATERIAL") == std::string::npos) {
				// 遇到其他关键字，表示材料块结束
				break;
			}
		}
		
		// 添加材料到inp_data
		inp_data.materials.push_back(material);
		
		// 同时将材料属性设置到模型中
		exdyna.mMatElastic.rho = material.density;
		exdyna.mMatElastic.E = material.E;
		exdyna.mMatElastic.v = material.v;
		exdyna.mMatElastic.update();  // 更新派生属性
		
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
		if (!currentPartPtr) return false;
		while (fin.peek() == '*') std::getline(fin, str);
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) continue;
			std::istringstream iss(str);
			std::string token;
			std::getline(iss, token, ','); // 单元ID
			std::vector<Types::Vertex_index> verticesIndex;
			while (std::getline(iss, token, ',')) {
				verticesIndex.push_back(std::stoul(token));
			}
			for (auto& idx : verticesIndex) idx -= 1;
			currentPartPtr->elements.push_back(verticesIndex);
		}
		return true;
	}

	/**
	 * @brief 处理节点集合数据
	 * 
	 * 解析输入文件中的节点集合定义，用于后续引用特定的节点组
	 * 支持处理实例引用，如"instance=Part-1-1"
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::SET_NODE_LIST() {
		std::string set_name, instance_name;
		parseSetWithInstance(str, set_name, instance_name);
		bool is_generate = (str.find("GENERATE") != std::string::npos);
		std::vector<std::size_t> set_nodes;
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) continue;
			std::istringstream iss(str);
			std::string token;
			while (!iss.eof()) {
				if (iss >> token) {
					if (token.back() == ',') token.pop_back();
					try { set_nodes.push_back(std::stoul(token)); } catch (...) {}
				}
				if (iss.peek() == ',') iss.ignore();
			}
		}
		if (is_generate && set_nodes.size() >= 2) {
			std::size_t start = set_nodes[0], end = set_nodes[1], step = (set_nodes.size() >= 3) ? set_nodes[2] : 1;
			set_nodes.clear();
			for (std::size_t i = start; i <= end; i += step) set_nodes.push_back(i);
		}
		for (auto& idx : set_nodes) idx -= 1;
		
		if (currentPartPtr && currentState == ParseState::PART_DEFINITION) {
			// 在Part定义中，直接添加到Part的节点集
			currentPartPtr->node_sets[set_name] = set_nodes;
		} 
		else if (currentAssemblyPtr && currentState == ParseState::ASSEMBLY_DEFINITION) {
			// 在Assembly定义中，需要考虑实例
			if (!instance_name.empty()) {
				// 如果指定了实例，添加到Assembly的节点集中特定实例下
				currentAssemblyPtr->node_sets[set_name][instance_name] = set_nodes;
				std::cout << "为实例 " << instance_name << " 添加节点集 " << set_name 
					      << " (包含 " << set_nodes.size() << " 个节点)" << std::endl;
			} 
			else {
				// 没有指定实例，视为装配级别的全局节点集 (旧行为保持兼容)
				std::map<std::string, std::vector<std::size_t>> empty_map;
				empty_map["__GLOBAL__"] = set_nodes; // 使用特殊键表示全局节点
				currentAssemblyPtr->node_sets[set_name] = empty_map;
				std::cout << "添加装配级全局节点集 " << set_name 
					      << " (包含 " << set_nodes.size() << " 个节点)" << std::endl;
			}
		}
		return true;
	}

	bool DataIn::BOUNDARY() {
		// 检查是否包含 op=NEW 参数，如果是，则清除当前步骤的相应边界条件
		bool op_new = false;
		if (str.find("OP=NEW") != std::string::npos) {
			op_new = true;
			
			// 如果在步骤定义中且有步骤数据，则清除当前步骤的相应边界条件并设置重置标志
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 根据边界条件类型（位移或速度）决定清除哪种边界条件
				if (str.find("TYPE=VELOCITY") != std::string::npos) {
					exdyna.steps.back().boundary.vel_nodes.clear();
					// 设置重置标志，表示本步骤重置了速度边界条件
					exdyna.steps.back().resetVelBoundary = true;
					std::cout << "清除步骤 " << exdyna.steps.back().name << " 的所有速度边界条件 (op=NEW)" << std::endl;
				} else {
					exdyna.steps.back().boundary.spc_nodes.clear();
					// 设置重置标志，表示本步骤重置了位移边界条件
					exdyna.steps.back().resetSpcBoundary = true;
					std::cout << "清除步骤 " << exdyna.steps.back().name << " 的所有位移边界条件 (op=NEW)" << std::endl;
				}
			} else {
				// 如果不在步骤定义中，则清除全局边界条件
				if (str.find("TYPE=VELOCITY") != std::string::npos) {
					exdyna.currentBoundary.vel_nodes.clear();
					std::cout << "清除所有全局速度边界条件 (op=NEW)" << std::endl;
				} else {
					exdyna.currentBoundary.spc_nodes.clear();
					std::cout << "清除所有全局位移边界条件 (op=NEW)" << std::endl;
				}
			}
		}

		// 处理不同类型的边界条件
		// 注意：当未指定type参数时（即 *Boundary 或 *Boundary, op=NEW 命令没有type参数），
		// 默认按SPC（位移）边界条件处理，调用SPC_NODE()方法
		
		// 检查当前命令行是否为单独的边界条件命令而没有后续数据
		// 这种情况在有连续的*Boundary命令时可能发生
		bool isEmptyBoundaryCommand = false;
		if (fin.peek() == '*') {
			// 下一行又是一个新的关键字，表示当前的*Boundary命令没有数据行
			isEmptyBoundaryCommand = true;
			
			// 即使是空命令，也需要记录已经处理过对应类型的重置操作
			if (op_new && currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 根据是否有velocity类型参数设置对应的重置标志
				if (str.find("TYPE=VELOCITY") != std::string::npos) {
					exdyna.steps.back().resetVelBoundary = true;
					std::cout << "遇到空的速度边界条件命令: " << str << " (仅执行重置操作)" << std::endl;
					std::cout << "为步骤 " << exdyna.steps.back().name << " 设置速度边界重置标志" << std::endl;
				} else {
					exdyna.steps.back().resetSpcBoundary = true;
					std::cout << "遇到空的位移边界条件命令: " << str << " (仅执行重置操作)" << std::endl;
					std::cout << "为步骤 " << exdyna.steps.back().name << " 设置位移边界重置标志" << std::endl;
				}
			} else {
				if (str.find("TYPE=VELOCITY") != std::string::npos) {
					std::cout << "遇到空的速度边界条件命令: " << str << " (未执行重置操作)" << std::endl;
				} else {
					std::cout << "遇到空的位移边界条件命令: " << str << " (未执行重置操作)" << std::endl;
				}
			}
		}
		
		// 只有当不是空的边界条件命令时，才调用相应的处理函数
		if (!isEmptyBoundaryCommand) {
			if (str.find("TYPE=VELOCITY") != std::string::npos) {
				VELOCITY_NODE();
			} else {
				SPC_NODE();
			}
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
		// 检查是否包含 op=NEW 参数，如果是，则清除当前步骤的相应载荷
		bool op_new = false;
		if (str.find("OP=NEW") != std::string::npos) {
			op_new = true;
			
			// 如果在步骤定义中且有步骤数据，则清除当前步骤的分布载荷并设置重置标志
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 清除分布载荷
				std::get<0>(exdyna.steps.back().gravity) = false;
				exdyna.steps.back().dsload.clear();
				
				// 设置重置标志
				exdyna.steps.back().resetDload = true;
				std::cout << "清除步骤 " << exdyna.steps.back().name << " 的所有分布载荷 (op=NEW)" << std::endl;
			}
		}
		
		// 从关键字行中提取幅值曲线名称
		std::string amp_name = extractNameFromKeyword(str, "AMPLITUDE=");
		
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
			Types::Real timePeriod = convertString<Types::Real>(token);
            
			// 如果当前在解析步骤，为当前步骤设置时间周期
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				exdyna.steps.back().timePeriod = timePeriod;
			}
            
			// 同时更新总时间
			exdyna.totalTime += timePeriod;
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

		// 如果下一行是另一个关键字（以*开头），表示这是一个纯粹的清除操作
		// 这种情况下，我们已经在BOUNDARY方法中处理了op=NEW参数，所以可以直接返回
		if (fin.peek() == '*') {
			std::cout << "位移边界条件已清除，无新添加的约束" << std::endl;
			return;
		}

		std::cout << "开始处理位移边界条件数据..." << std::endl;
		
		// 读取单点约束数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
            std::cout << "  解析位移边界行: " << str << std::endl;
			
			// 创建临时约束数据
			std::pair<std::vector<std::size_t>, std::array<std::size_t, 3>> boundary_data;
			std::istringstream iss(str);
			std::string token;

			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				// 使用预定义的节点集合
				std::getline(iss, token, ',');
				std::cout << "  使用节点集: " << token << std::endl;
				
				// 检查节点集是否存在
				if (exdyna.map_set_node_list.find(token) != exdyna.map_set_node_list.end()) {
					boundary_data.first = exdyna.map_set_node_list[token];
				} else {
					std::cerr << "  警告: 找不到节点集 " << token << std::endl;
					continue; // 跳过这一行
				}
			}
			else {
				// 单个节点
				std::getline(iss, token, ',');
				std::cout << "  使用单个节点: " << token << std::endl;
				try {
					boundary_data.first.push_back(std::stoi(token) - 1);
				} catch (const std::exception& e) {
					std::cerr << "  警告: 解析节点索引出错: " << e.what() << std::endl;
					continue; // 跳过这一行
				}
			}

			// 处理特殊约束类型
			if (str.find("PINNED") != std::string::npos || str.find("ENCASTRE") != std::string::npos) {
				std::cout << "  使用特殊约束类型: " << (str.find("PINNED") != std::string::npos ? "PINNED" : "ENCASTRE") << std::endl;
				boundary_data.second.at(0) = 0;  // 起始自由度
				boundary_data.second.at(1) = 2;  // 结束自由度
				boundary_data.second.at(2) = 1;  // 增量
			}
			else if (str.find("XSYMM") != std::string::npos) {
				// X对称：约束X方向位移
				std::cout << "  使用特殊约束类型: XSYMM (X对称)" << std::endl;
				boundary_data.second.at(0) = 0;  // X方向 (第1个自由度)
				boundary_data.second.at(1) = 0;  // 仅X方向
				boundary_data.second.at(2) = 1;  // 增量
			}
			else if (str.find("YSYMM") != std::string::npos) {
				// Y对称：约束Y方向位移
				std::cout << "  使用特殊约束类型: YSYMM (Y对称)" << std::endl;
				boundary_data.second.at(0) = 1;  // Y方向 (第2个自由度)
				boundary_data.second.at(1) = 1;  // 仅Y方向
				boundary_data.second.at(2) = 1;  // 增量
			}
			else if (str.find("ZSYMM") != std::string::npos) {
				// Z对称：约束Z方向位移
				std::cout << "  使用特殊约束类型: ZSYMM (Z对称)" << std::endl;
				boundary_data.second.at(0) = 2;  // Z方向 (第3个自由度)
				boundary_data.second.at(1) = 2;  // 仅Z方向
				boundary_data.second.at(2) = 1;  // 增量
			}
			else {
				// 读取自由度约束定义
				std::getline(iss, token, ',');
				if (token.empty()) {
					// 如果没有指定自由度，默认约束所有自由度
					std::cout << "  未指定自由度，默认约束所有自由度 (1-3)" << std::endl;
					boundary_data.second.at(0) = 0;
					boundary_data.second.at(1) = 2;
					boundary_data.second.at(2) = 1;
				} else {
					std::cout << "  自由度起始值: " << token << std::endl;
					boundary_data.second.at(0) = std::stoi(token) - 1;  // 起始自由度
					
					std::getline(iss, token, ',');
					if (token.empty()) {
						// 如果没有指定结束自由度，则设置为与起始自由度相同
						boundary_data.second.at(1) = boundary_data.second.at(0);
						std::cout << "  自由度结束值: 与起始值相同" << std::endl;
					}
					else {
						boundary_data.second.at(1) = std::stoi(token) - 1;  // 结束自由度
						std::cout << "  自由度结束值: " << token << std::endl;
					}
					
					// 读取增量
					std::getline(iss, token, ',');
					if (token.empty()) {
						boundary_data.second.at(2) = 1;  // 默认增量为1
						std::cout << "  增量: 默认为1" << std::endl;
					} else {
						boundary_data.second.at(2) = std::stoi(token);
						std::cout << "  增量: " << token << std::endl;
					}
				}
			}
            
			// 添加到当前步骤的约束列表中
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 只添加到当前解析中的步骤
				exdyna.steps.back().boundary.spc_nodes.push_back(boundary_data);
				std::cout << "在步骤 " << exdyna.steps.back().name << " 添加固定约束 (共 " 
                          << boundary_data.first.size() << " 个节点)" << std::endl;
			} else {
				// 如果不在Step定义中，则添加到全局约束
				exdyna.currentBoundary.spc_nodes.push_back(boundary_data);
				std::cout << "添加全局固定约束 (共 " << boundary_data.first.size() << " 个节点)" << std::endl;
			}
		}
		
		std::cout << "位移边界条件处理完成" << std::endl;
	}

	/**
	 * @brief 处理初始速度数据
	 * 
	 * 解析输入文件中的初始速度数据，为节点设置初始速度
	 * 支持处理实例引用和节点集合
	 */
	void DataIn::INITIAL_VELOCITY() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*')
			std::getline(fin, str);
		
		// 读取初始速度数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) continue;
			
			toUpperCase(str);
			exdyna.ini_vel_generation.emplace_back();
			std::istringstream iss(str);
			std::string token;

			// 处理节点引用（节点集合或单个节点）
			std::string node_ref;
			std::getline(iss, token, ',');
			node_ref = token;
			removeSpaces(node_ref);
			
			// 检查是否是节点集引用
			bool is_set = (node_ref.find("SET-") != std::string::npos || node_ref.find("NSET-") != std::string::npos);
			
			// 处理可能的实例前缀，如"Part-1-1.Set-1"
			std::string instance_name = "";
			size_t dot_pos = node_ref.find('.');
			if (dot_pos != std::string::npos) {
				instance_name = node_ref.substr(0, dot_pos);
				node_ref = node_ref.substr(dot_pos + 1);
				std::cout << "引用实例 " << instance_name << " 中的";
			}
			
			// 处理节点集合或单个节点
			if (is_set) {
				// 使用预定义的节点集合
				std::string set_name = node_ref;
				// 构建完整的集合名称（包含实例前缀，如果有）
				if (!instance_name.empty()) {
					set_name = instance_name + "." + node_ref;
				}
				
				// 检查节点集是否存在
				if (exdyna.map_set_node_list.find(set_name) != exdyna.map_set_node_list.end()) {
					exdyna.ini_vel_generation.back().first = exdyna.map_set_node_list[set_name];
					std::cout << "为节点集 " << set_name << " 设置初始速度" << std::endl;
				} else {
					std::cerr << "警告: 找不到节点集 " << set_name << "，跳过初始速度设置" << std::endl;
					exdyna.ini_vel_generation.pop_back();
					continue;
				}
			}
			else {
				// 单个节点，转换索引（从1到0）
				try {
					int node_index = std::stoi(node_ref) - 1;
					exdyna.ini_vel_generation.back().first.push_back(node_index);
					std::cout << "为节点 " << (node_index + 1) << " 设置初始速度" << std::endl;
				}
				catch (const std::exception& e) {
					std::cerr << "警告: 解析节点索引时出错: " << node_ref << std::endl;
					exdyna.ini_vel_generation.pop_back();
					continue;
				}
			}

			// 读取自由度索引和速度值
			try {
				// 自由度索引
				if (std::getline(iss, token, ',')) {
					exdyna.ini_vel_generation.back().second.first = std::stoi(token) - 1;  // 自由度索引（从1到0）
				} else {
					// 默认为第一个自由度
					exdyna.ini_vel_generation.back().second.first = 0;
				}
				
				// 速度值
				if (std::getline(iss, token, ',')) {
					exdyna.ini_vel_generation.back().second.second = convertString<Types::Real>(token);
					std::cout << "  - 自由度: " << (exdyna.ini_vel_generation.back().second.first + 1) 
					          << ", 速度: " << exdyna.ini_vel_generation.back().second.second << std::endl;
				} else {
					// 如果没有指定速度值，默认为0
					exdyna.ini_vel_generation.back().second.second = 0.0;
					std::cout << "  - 自由度: " << (exdyna.ini_vel_generation.back().second.first + 1) 
					          << ", 速度: 0.0 (默认)" << std::endl;
				}
			}
			catch (const std::exception& e) {
				std::cerr << "警告: 解析速度数据时出错: " << e.what() << std::endl;
				exdyna.ini_vel_generation.pop_back();
			}
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
		
		// 如果下一行是另一个关键字（以*开头），表示这是一个纯粹的清除操作
		// 这种情况下，我们已经在BOUNDARY方法中处理了op=NEW参数，所以可以直接返回
		if (fin.peek() == '*') {
			std::cout << "速度边界条件已清除，无新添加的约束" << std::endl;
			return;
		}
		
		std::cout << "开始处理速度边界条件数据..." << std::endl;
		
		// 读取速度边界条件数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			std::cout << "  解析速度边界行: " << str << std::endl;
			
			std::istringstream iss(str);
			std::string token;

			// 创建临时结构存储速度边界条件
			std::pair<std::array<std::size_t, 2>, Types::Real> temp_pair;
			temp_pair.first[0] = 0;
			temp_pair.first[1] = 0;
			temp_pair.second = 0.0;
			
			std::pair<std::vector<std::size_t>, std::pair<std::array<std::size_t, 2>, Types::Real>> boundary_data;
			boundary_data.second = temp_pair;
			
			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				std::getline(iss, token, ',');
				std::cout << "  使用节点集: " << token << std::endl;
				
				// 检查节点集是否存在
				if (exdyna.map_set_node_list.find(token) != exdyna.map_set_node_list.end()) {
					boundary_data.first = exdyna.map_set_node_list[token];
				} else {
					std::cerr << "  警告: 找不到节点集 " << token << std::endl;
					continue; // 跳过这一行
				}
			}
			else {
				std::getline(iss, token, ',');
				std::cout << "  使用单个节点: " << token << std::endl;
				try {
					std::vector<std::size_t> temp_vector;
					temp_vector.emplace_back(std::stoi(token) - 1);
					boundary_data.first = temp_vector;
				} catch (const std::exception& e) {
					std::cerr << "  警告: 解析节点索引出错: " << e.what() << std::endl;
					continue; // 跳过这一行
				}
			}

			// 读取自由度信息
			std::getline(iss, token, ',');
			if (token.empty()) {
				std::cerr << "  警告: 缺少自由度信息" << std::endl;
				continue; // 跳过这一行
			}
			
			std::cout << "  自由度起始值: " << token << std::endl;
			boundary_data.second.first[0] = std::stoi(token) - 1;
			
			std::getline(iss, token, ',');
			if (token.empty()) {
				boundary_data.second.first[1] = boundary_data.second.first[0];
				std::cout << "  自由度结束值: 与起始值相同" << std::endl;
			}
			else {
				boundary_data.second.first[1] = std::stoi(token) - 1;
				std::cout << "  自由度结束值: " << token << std::endl;
			}

			// 读取速度值
			if (std::getline(iss, token, ',')) {
				boundary_data.second.second = convertString<Types::Real>(token);
				std::cout << "  速度值: " << boundary_data.second.second << std::endl;
			} else {
				std::cout << "  默认速度值: 0.0" << std::endl;
			}
            
			// 添加到当前步骤的约束列表中
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 只添加到当前解析中的步骤
				exdyna.steps.back().boundary.vel_nodes.push_back(boundary_data);
				std::cout << "在步骤 " << exdyna.steps.back().name << " 添加速度约束 (共 " 
                          << boundary_data.first.size() << " 个节点)" << std::endl;
			} else {
				// 如果不在Step定义中，则添加到全局约束
				exdyna.currentBoundary.vel_nodes.push_back(boundary_data);
				std::cout << "添加全局速度约束 (共 " << boundary_data.first.size() << " 个节点)" << std::endl;
			}
		}
		
		std::cout << "速度边界条件处理完成" << std::endl;
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

		// 设置重力参数，根据当前解析状态决定存储位置
		if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
			// 存储在当前步骤的gravity中
			std::get<0>(exdyna.steps.back().gravity) = true;
			std::get<1>(exdyna.steps.back().gravity) = amp_name;
			std::get<2>(exdyna.steps.back().gravity) = value;
			std::get<3>(exdyna.steps.back().gravity) = direction_x;
			std::get<4>(exdyna.steps.back().gravity) = direction_y;
			std::get<5>(exdyna.steps.back().gravity) = direction_z;
			std::cout << "为步骤 " << exdyna.steps.back().name << " 设置重力加速度: " 
					  << value << " 方向: (" << direction_x << ", " << direction_y << ", " << direction_z << ")" << std::endl;
		} else {
			// 否则存储在全局变量中（不应该发生，但作为后备）
			std::cout << "警告：在步骤外定义重力，这可能导致不一致的行为" << std::endl;
		}
	}

	/**
	 * @brief 处理幅值曲线数据
	 * 
	 * 解析输入文件中的幅值曲线定义，用于后续引用特定的时间幅值曲线
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::AMP() {
		// 从关键字行中提取幅值曲线名称
		std::string amp_name = extractNameFromKeyword(str);

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

	/**
	 * @brief 处理单元集合数据
	 * 
	 * 解析输入文件中的单元集合定义，用于后续引用特定的单元组
	 * 支持处理实例引用，如"instance=Part-1-1"
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::SET_ELE_LIST() {
		std::string set_name, instance_name;
		parseSetWithInstance(str, set_name, instance_name);
		std::vector<std::size_t> element_indices;
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) continue;
			std::istringstream iss(str);
			std::string token;
			std::vector<std::size_t> line_values;
			while (!iss.eof()) {
				if (iss >> token) {
					if (token.back() == ',') token.pop_back();
					try { line_values.push_back(std::stoul(token)); } catch (...) {}
				}
				if (iss.peek() == ',') iss.ignore();
			}
			if (line_values.size() == 3) {
				std::size_t start = line_values[0], end = line_values[1], step = line_values[2];
				for (std::size_t i = start; i <= end; i += step) element_indices.push_back(i);
			} else {
				element_indices.insert(element_indices.end(), line_values.begin(), line_values.end());
			}
		}
		for (auto& idx : element_indices) idx -= 1;
		
		if (currentPartPtr && currentState == ParseState::PART_DEFINITION) {
			// 在Part定义中，直接添加到Part的单元集
			currentPartPtr->element_sets[set_name] = element_indices;
		} 
		else if (currentAssemblyPtr && currentState == ParseState::ASSEMBLY_DEFINITION) {
			// 在Assembly定义中，需要考虑实例
			if (!instance_name.empty()) {
				// 如果指定了实例，添加到Assembly的单元集中特定实例下
				currentAssemblyPtr->element_sets[set_name][instance_name] = element_indices;
				std::cout << "为实例 " << instance_name << " 添加单元集 " << set_name 
					      << " (包含 " << element_indices.size() << " 个单元)" << std::endl;
			} 
			else {
				// 没有指定实例，视为装配级别的全局单元集 (旧行为保持兼容)
				std::map<std::string, std::vector<std::size_t>> empty_map;
				empty_map["__GLOBAL__"] = element_indices; // 使用特殊键表示全局单元
				currentAssemblyPtr->element_sets[set_name] = empty_map;
				std::cout << "添加装配级全局单元集 " << set_name 
					      << " (包含 " << element_indices.size() << " 个单元)" << std::endl;
			}
		}
		return true;
	}

	/**
	 * @brief 处理表面数据
	 * 
	 * 解析输入文件中的表面定义，用于后续引用特定的表面
	 * 支持处理实例引用
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::SURFACE() {
		// 从关键字行中提取表面名称和实例名称（如果有）
		std::string surf_name = extractNameFromKeyword(str, "NAME=");
		std::string instance_name = extractNameFromKeyword(str, "INSTANCE=");
		
		// 处理TYPE参数，如果有的话
		std::string surf_type_global = extractNameFromKeyword(str, "TYPE=");
		
		if (surf_name.empty()) {
			std::cerr << "警告: 无法解析表面名称" << std::endl;
			// 跳过这个表面定义的数据
			while (fin.peek() != '*' && fin.peek() != EOF) {
				std::getline(fin, str);
			}
			return false;
		}
		
		// 标准化表面名称
		removeSpaces(surf_name);
		
		// 如果有实例名称，加上前缀
		std::string full_surf_name = surf_name;
		if (!instance_name.empty()) {
			full_surf_name = instance_name + "." + surf_name;
			std::cout << "在实例 " << instance_name << " 中定义表面: " << surf_name << std::endl;
		} else {
			std::cout << "定义表面: " << surf_name << std::endl;
		}
		
		// 逐行读取表面定义
		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);
			if (str.empty()) continue;
			
			toUpperCase(str);
			std::istringstream iss(str);
			std::string token;
			
			// 读取单元集合引用
			std::string ele_set_name;
			std::getline(iss, token, ',');
			ele_set_name = token;
			removeSpaces(ele_set_name);
			
			// 处理可能的实例前缀，如"Part-1-1.Set-1"
			std::string element_instance = instance_name;
			size_t dot_pos = ele_set_name.find('.');
			if (dot_pos != std::string::npos) {
				element_instance = ele_set_name.substr(0, dot_pos);
				ele_set_name = ele_set_name.substr(dot_pos + 1);
			}
			
			// 读取表面类型
			std::string surf_face_type;
			std::getline(iss, token, ',');
			surf_face_type = token;
			removeSpaces(surf_face_type);
			
			// 应用全局表面类型，如果没有指定局部类型的话
			if (surf_face_type.empty() && !surf_type_global.empty()) {
				surf_face_type = surf_type_global;
			}
			
			// 使用实例名称构建完整的单元集合名称
			std::string full_ele_set_name = ele_set_name;
			if (!element_instance.empty()) {
				full_ele_set_name = element_instance + "." + ele_set_name;
			}
			
			// 添加到表面映射
			std::pair<std::string, std::string> surf_definition(full_ele_set_name, surf_face_type);
			exdyna.map_set_surface_list[full_surf_name] = surf_definition;
			
			std::cout << "  - 表面组成: 单元集=" << full_ele_set_name << ", 面类型=" << surf_face_type << std::endl;
		}
		
		return true;
	}

	/**
	 * @brief 处理分布面载荷数据
	 * 
	 * 解析输入文件中的分布面载荷数据
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::DSLOAD() {
		// 检查是否包含 op=NEW 参数，如果是，则清除当前步骤的相应载荷
		bool op_new = false;
		if (str.find("OP=NEW") != std::string::npos) {
			op_new = true;
			
			// 如果在步骤定义中且有步骤数据，则清除当前步骤的分布面载荷并设置重置标志
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 清除分布面载荷
				exdyna.steps.back().dsload.clear();
				
				// 设置重置标志
				exdyna.steps.back().resetDsload = true;
				std::cout << "清除步骤 " << exdyna.steps.back().name << " 的所有分布面载荷 (op=NEW)" << std::endl;
			}
		}

		// 从关键字行中提取幅值曲线名称
		std::string amp_name = extractNameFromKeyword(str, "AMPLITUDE=");

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
			
			// 创建新的分布面载荷元组
			std::tuple<std::string, std::string, std::string, Types::Real> dsload_entry;
			std::get<0>(dsload_entry) = surf_name;
			std::get<1>(dsload_entry) = load_type;
			std::get<2>(dsload_entry) = amp_name;
			std::get<3>(dsload_entry) = load_value;
			
			// 根据当前解析状态决定存储位置
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 存储在当前步骤的dsload中
				exdyna.steps.back().dsload.push_back(dsload_entry);
				std::cout << "为步骤 " << exdyna.steps.back().name 
						  << " 添加分布面载荷: 表面=" << surf_name 
						  << " 类型=" << load_type 
						  << " 值=" << load_value << std::endl;
			} else {
				// 否则存储在全局变量中（不应该发生，但作为后备）
				std::cout << "警告：在步骤外定义分布面载荷，这可能导致不一致的行为" << std::endl;
			}
		}
		return true;
	}

	bool DataIn::OUTPUT_INTERVAL() {
		// 跳过可能存在的注释行
		while (fin.peek() == '*')
			std::getline(fin, str);
		
		// 读取输出间隔值
		if (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			std::istringstream iss(str);
			std::string token;
			
			// 读取输出间隔
			std::getline(iss, token, ',');
			if (!token.empty()) {
				// 将读取的值设置为exdyna的time_interval
				exdyna.time_interval = convertString<Types::Real>(token);
				std::cout << "从inp文件中读取到输出间隔: " << exdyna.time_interval << std::endl;
				
				// 标记time_interval已被设置
				exdyna.time_interval_set = true;
			}
		}
		
		return true;
	}

	bool DataIn::OUTPUT() {
		// 读取Abaqus的*Output关键字行
		std::string original_line = str;
		
		// 转换为大写以便查找参数（不区分大小写）
		std::string upper_line = original_line;
		toUpperCase(upper_line);
		
		// 检查是否包含number interval参数
		size_t pos = upper_line.find("NUMBER INTERVAL=");
		if (pos != std::string::npos) {
			// 提取number interval的值
			std::string interval_str = original_line.substr(pos + 16); // 16是"NUMBER INTERVAL="的长度
			// 查找逗号或空格，确定值的结束位置
			size_t end_pos = interval_str.find_first_of(", ");
			if (end_pos != std::string::npos) {
				interval_str = interval_str.substr(0, end_pos);
			}
			
			// 转换为数值
			try {
				int number_interval = std::stoi(interval_str);
				
				// 根据Abaqus的number interval设置time_interval
				// number interval是总步数的划分，所以time_interval = totalTime / number_interval
				if (number_interval > 0) {
					exdyna.time_interval = exdyna.totalTime / number_interval;
					exdyna.time_interval_set = true;
					std::cout << "从inp文件的*Output设置中读取到number interval=" << number_interval 
					          << "，设置time_interval=" << exdyna.time_interval << std::endl;
				}
			} catch (const std::exception& e) {
				std::cerr << "解析OUTPUT的number interval参数时出错: " << e.what() << std::endl;
			}
		}
		
		// 跳过*Output下面的内容，直到下一个关键字或文件结束
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
		}
		
		return true;
	}

	// 新增：transferToExDyna方法的实现
	void DataIn::transferToExDyna() {
		std::cout << "\n=== 开始转换数据到exDyna3D结构 ===" << std::endl;
		
		// 统计数量
		size_t totalNodes = 0;
		size_t totalElements = 0;
		
		// 1. 统计所有Part的节点和单元数量
		std::cout << "处理的部件数量: " << inp_data.parts.size() << std::endl;
		for (const auto& part : inp_data.parts) {
			totalNodes += part.nodes.size();
			totalElements += part.elements.size();
			std::cout << "- 部件 '" << part.name << "': " 
					  << part.nodes.size() << " 个节点, " 
					  << part.elements.size() << " 个单元" << std::endl;
		}
		
		// 如果没有Part数据（即inp_data未被填充），则直接返回，不进行转换
		if (inp_data.parts.empty()) {
			std::cout << "警告：没有找到有效的Part数据，跳过转换。" << std::endl;
			return;
		}
		
		// 2. 为exDyna3D分配空间
		exdyna.vertices.clear();
		exdyna.hexahedron_elements.clear();
		exdyna.vertices.reserve(totalNodes);
		exdyna.hexahedron_elements.reserve(totalElements);
		
		// 3. 首先处理所有Part的数据
		std::map<std::string, size_t> partNodeOffset; // 记录每个Part节点的全局偏移量
		std::map<std::string, size_t> partElementOffset; // 记录每个Part单元的全局偏移量
		size_t globalNodeIndex = 0;
		size_t globalElementIndex = 0;
		
		// 3.1 添加节点到exdyna
		for (const auto& part : inp_data.parts) {
			// 记录当前Part的节点偏移量
			partNodeOffset[part.name] = globalNodeIndex;
			
			// 添加节点到exdyna
			for (const auto& node : part.nodes) {
				exdyna.vertices.push_back(node);
				globalNodeIndex++;
			}
		}
		
		// 3.2 添加单元到exdyna
		for (const auto& part : inp_data.parts) {
			// 记录当前Part的单元偏移量
			partElementOffset[part.name] = globalElementIndex;
			size_t partOffset = partNodeOffset[part.name];
			
			// 添加单元到exdyna
			for (const auto& element : part.elements) {
				exdyna.hexahedron_elements.emplace_back();
				auto& newElement = exdyna.hexahedron_elements.back();
				
				// 设置单元的物理ID和材料ID（默认为1）
				newElement.set_PID(1);
				newElement.set_MID(1);
				
				// 复制单元的节点索引，并调整为全局索引
				std::vector<Types::Vertex_index> globalNodes;
				globalNodes.reserve(element.size());
				
				for (const auto& localIdx : element) {
					globalNodes.push_back(localIdx + partOffset);
				}
				
				// 设置调整后的节点索引
				newElement.set_verticesIndex(globalNodes);
				globalElementIndex++;
			}
		}
		
		// 4. 处理装配中的实例，更新实例的索引信息
		if (!inp_data.assemblies.empty()) {
			std::cout << "处理的装配数量: " << inp_data.assemblies.size() << std::endl;
			
			for (auto& assembly : inp_data.assemblies) {
				std::cout << "- 处理装配 '" << assembly.name << "' 及其实例" << std::endl;
				
				// 处理实例
				for (auto& instance : assembly.instances) {
					std::cout << "  - 实例 '" << instance.name << "' 引用部件 '" << instance.part_name << "'" << std::endl;
					
					// 找到对应的部件
					auto partIt = std::find_if(inp_data.parts.begin(), inp_data.parts.end(), 
						[&](const InpPart& p) { return p.name == instance.part_name; });
					
					if (partIt != inp_data.parts.end()) {
						// 设置实例的节点和单元起始索引
						instance.node_start_index = partNodeOffset[instance.part_name];
						instance.element_start_index = partElementOffset[instance.part_name];
						
						std::cout << "    节点起始索引: " << instance.node_start_index 
						          << ", 单元起始索引: " << instance.element_start_index << std::endl;
						
						// 处理坐标变换（平移）
						// 注意：这里简化处理，假设不需要坐标变换，在实际应用中可以添加
					} else {
						std::cerr << "  警告: 找不到部件 '" << instance.part_name << "'" << std::endl;
					}
				}
				
				// 5. 处理装配级别的节点集和单元集
				// 5.1 处理节点集
				for (const auto& nodeSetPair : assembly.node_sets) {
					const std::string& setName = nodeSetPair.first;
					const auto& instanceMap = nodeSetPair.second;
					
					// 合并所有实例的节点集到一个全局节点集
					std::vector<std::size_t> globalNodeSet;
					
					for (const auto& instPair : instanceMap) {
						const std::string& instName = instPair.first;
						const auto& localIndices = instPair.second;
						
						if (instName == "__GLOBAL__") {
							// 全局节点索引直接添加
							globalNodeSet.insert(globalNodeSet.end(), localIndices.begin(), localIndices.end());
						} else {
							// 查找实例
							auto instIt = std::find_if(assembly.instances.begin(), assembly.instances.end(),
								[&](const InpInstance& inst) { return inst.name == instName; });
								
							if (instIt != assembly.instances.end()) {
								// 添加转换后的节点索引
								for (const auto& localIdx : localIndices) {
									globalNodeSet.push_back(localIdx + instIt->node_start_index);
								}
							}
						}
					}
					
					// 添加到exdyna的节点集映射
					exdyna.map_set_node_list[setName] = globalNodeSet;
					std::cout << "  添加装配级节点集: " << setName << " 包含 " << globalNodeSet.size() << " 个节点" << std::endl;
				}
				
				// 5.2 处理单元集
				for (const auto& elemSetPair : assembly.element_sets) {
					const std::string& setName = elemSetPair.first;
					const auto& instanceMap = elemSetPair.second;
					
					// 合并所有实例的单元集到一个全局单元集
					std::vector<std::size_t> globalElemSet;
					
					for (const auto& instPair : instanceMap) {
						const std::string& instName = instPair.first;
						const auto& localIndices = instPair.second;
						
						if (instName == "__GLOBAL__") {
							// 全局单元索引直接添加
							globalElemSet.insert(globalElemSet.end(), localIndices.begin(), localIndices.end());
						} else {
							// 查找实例
							auto instIt = std::find_if(assembly.instances.begin(), assembly.instances.end(),
								[&](const InpInstance& inst) { return inst.name == instName; });
								
							if (instIt != assembly.instances.end()) {
								// 添加转换后的单元索引
								for (const auto& localIdx : localIndices) {
									globalElemSet.push_back(localIdx + instIt->element_start_index);
								}
							}
						}
					}
					
					// 添加到exdyna的单元集映射
					exdyna.map_set_ele_list[setName] = globalElemSet;
					std::cout << "  添加装配级单元集: " << setName << " 包含 " << globalElemSet.size() << " 个单元" << std::endl;
				}
			}
		}
		
		// 6. 处理材料属性
		if (!inp_data.materials.empty()) {
			std::cout << "处理材料数据: " << inp_data.materials.size() << " 个材料" << std::endl;
			
			// 这里只取第一个材料（假设单一材料）
			if (!inp_data.materials.empty()) {
				const auto& material = inp_data.materials[0];
				exdyna.mMatElastic.rho = material.density;
				exdyna.mMatElastic.E = material.E;
				exdyna.mMatElastic.v = material.v;
				std::cout << "设置材料属性: 密度=" << material.density 
						  << " E=" << material.E 
						  << " v=" << material.v << std::endl;
			}
		}
		
		std::cout << "数据转换完成: " 
				  << exdyna.vertices.size() << " 个节点, " 
				  << exdyna.hexahedron_elements.size() << " 个单元" 
				  << std::endl;
		std::cout << "=== 转换完成 ===\n" << std::endl;
	}
}

