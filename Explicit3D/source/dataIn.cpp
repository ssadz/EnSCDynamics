#include"../include/dataIn.h"
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include <string>
#include <cctype>
#include "../include/InpData.h"
#include "spdlog/spdlog.h"


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
		
		// 检查下一行是否包含平移信息
		std::streampos currentPos = fin.tellg(); // 保存当前位置
		if (fin.peek() != '*' && !fin.eof()) {
			std::string nextLine;
			std::getline(fin, nextLine);
			
			// 检查下一行是否包含可能的平移量数据（三个浮点数）
			std::istringstream iss(nextLine);
			Types::Real x, y, z;
			bool isTranslation = false;
			
			// 尝试解析三个浮点数
			char peek_char = iss.peek();
			while (std::isspace(peek_char) && peek_char != EOF) {
				iss.ignore(1);
				peek_char = iss.peek();
			}
			
			if (iss >> x) {
				// 检查第一个数字后是否有逗号，如果有则移除
				if (iss.peek() == ',') iss.ignore(1);
				if (iss >> y) {
					// 检查第二个数字后是否有逗号，如果有则移除
					if (iss.peek() == ',') iss.ignore(1);
					if (iss >> z) {
						isTranslation = true;
					}
				}
			}
			
			if (isTranslation) {
				// 设置平移量
				inst.translation[0] = x;
				inst.translation[1] = y;
				inst.translation[2] = z;
				spdlog::debug("  实例 {} 的平移量: ({}, {}, {})", instance_name, x, y, z);
			} else {
				// 不是平移量数据，回到之前的位置
				fin.seekg(currentPos);
			}
		}
		
		currentAssemblyPtr->instances.push_back(inst);
		currentInstanceName = instance_name;
		currentInstancePart = part_name;
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
		spdlog::debug("结束部件定义: {}", currentPartName); // 替换 std::cout
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
		spdlog::debug("结束装配定义: {}", currentAssemblyName); // 替换 std::cout
		return true;
	}

	/**
	 * @brief 处理步骤结束标记
	 * 
	 * @return 处理成功返回true
	 */
	bool DataIn::END_STEP() {
		currentState = ParseState::GLOBAL;
		spdlog::debug("结束步骤定义: {}", currentStepName); // 替换 std::cout
		
		// 增加更多调试信息
		if (!exdyna.steps.empty()) {
			StepData& currentStep = exdyna.steps.back();
			spdlog::debug("------------------步骤调试信息------------------"); // 替换 std::cout
			spdlog::debug("步骤名称: {}", currentStep.name); // 替换 std::cout
			spdlog::debug("时间周期: {}", currentStep.timePeriod); // 替换 std::cout
			spdlog::debug("本步骤定义的固定约束节点条件数: {}", currentStep.boundary.spc_nodes.size()); // 替换 std::cout
			spdlog::debug("本步骤定义的速度约束节点条件数: {}", currentStep.boundary.vel_nodes.size()); // 替换 std::cout
			
			// 打印重置标记信息
			if (currentStep.resetSpcBoundary) {
				spdlog::debug("本步骤重置了所有位移边界条件 (Boundary, op=NEW)"); // 替换 std::cout
			}
			if (currentStep.resetVelBoundary) {
				spdlog::debug("本步骤重置了所有速度边界条件 (Boundary, op=NEW, type=VELOCITY)"); // 替换 std::cout
			}
			if (currentStep.resetDload) {
				spdlog::debug("本步骤重置了所有分布载荷 (DLOAD, op=NEW)"); // 替换 std::cout
			}
			if (currentStep.resetDsload) {
				spdlog::debug("本步骤重置了所有分布面载荷 (DSLOAD, op=NEW)"); // 替换 std::cout
			}
			
			// 提示关于边界条件继承的情况
			if (exdyna.steps.size() > 1) {
				if (!currentStep.resetSpcBoundary && !currentStep.resetVelBoundary &&
					!currentStep.resetDload && !currentStep.resetDsload) {
					spdlog::debug("注意: 在运行时，如果本步骤未定义某类条件，将从前一步骤继承"); // 替换 std::cout
				} else {
					spdlog::debug("注意: 本步骤有重置标记，部分条件将不会从前面步骤继承"); // 替换 std::cout
				}
			}
			
			spdlog::debug("------------------------------------------------"); // 替换 std::cout
		} else {
			spdlog::warn("当前步骤数据未保存!"); // 替换 std::cout
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
		spdlog::debug("结束实例定义"); // 替换 std::cout
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
			spdlog::error("无法打开文件！"); // 替换 std::cout
			exit(0);
		}
		else {
			spdlog::info("成功打开文件 {}！", fileName); // 替换 std::cout
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

		// 创建初始步骤"Initial"，用于收集*Step外的边界条件
		StepData initialStep;
		initialStep.name = "Initial";
		initialStep.timePeriod = 0.0; // 初始步骤的计算时间为0
		exdyna.steps.push_back(initialStep);
		spdlog::debug("创建初始步骤 'Initial' 用于收集*Step外的条件"); // 替换 std::cout

		// 逐行读取并处理文件内容
		while (!fin.eof()) {
			if (fin.peek() == '*') {  // 检查是否是关键字行
				std::getline(fin, str);
				if (str.length() > 1) {
					char secondChar = str[1];
					if (secondChar != '*') {  // 忽略注释行（以 ** 开头）
						toUpperCase(str);  // 转换为大写，方便匹配
						spdlog::debug("处理关键字行: {}", str);  // 添加调试输出, 替换 std::cout
						
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
							spdlog::debug("进入部件定义模式: {}", currentPartName); // 替换 std::cout
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
							spdlog::debug("进入装配定义模式: {}", currentAssemblyName); // 替换 std::cout
							continue;
						}
						else if (str.find("*MATERIAL") != std::string::npos) {
							currentMaterialName = extractNameFromKeyword(str);
							currentState = ParseState::MATERIAL_DEFINITION;
							spdlog::debug("进入材料定义模式: {}", currentMaterialName); // 替换 std::cout
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
							spdlog::debug("进入步骤定义模式: {}", currentStepName); // 替换 std::cout
							continue;
						}
						else if (str.find("*INITIAL CONDITIONS") != std::string::npos) {
							currentState = ParseState::PREDEFINED_FIELD;
							spdlog::debug("进入预定义场模式"); // 替换 std::cout
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
									spdlog::debug("为单元集 {} 应用材料 {}", elset, material); // 替换 std::cout
									
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
										spdlog::debug("设置材料 {} 的密度: {}", currentMaterialName, density); // 替换 std::cout
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
										
										spdlog::debug("设置材料 {} 的弹性模量: {}，泊松比: {}", currentMaterialName, elasticModulus, poissonRatio); // 替换 std::cout
										
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
									spdlog::debug("**特殊处理步骤内边界条件命令: {}", str); // 替换 std::cout
									BOUNDARY();  // 直接调用处理函数
									spdlog::debug("**步骤内边界条件命令处理完成"); // 替换 std::cout
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
											spdlog::debug("从输出设置获取间隔数: {}, 设置time_interval={}", interval, exdyna.time_interval); // 替换 std::cout
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
		spdlog::debug("成功读取文件 {}！", fileName); // 替换 std::cout
		
		// 新增：调用transferToExDyna()将inp_data转换为exdyna数据
		transferToExDyna();
		
		// 更新材料属性，确保在所有数据读取后执行
		exdyna.mMatElastic.update();
		
		// 打印解析统计信息
		spdlog::debug("解析完成: "); // 替换 std::cout
		if (!currentPartName.empty()) {
			spdlog::debug("- 部件名称: {}", currentPartName); // 替换 std::cout
		}
		if (!currentAssemblyName.empty()) {
			spdlog::debug("- 装配名称: {}", currentAssemblyName); // 替换 std::cout
		}
		if (!currentStepName.empty()) {
			spdlog::debug("- 步骤名称: {}", currentStepName); // 替换 std::cout
		}
		if (!currentMaterialName.empty()) {
			spdlog::debug("- 材料名称: {}", currentMaterialName); // 替换 std::cout
		}
		spdlog::debug("- 节点数量: {}", exdyna.vertices.size()); // 替换 std::cout
		spdlog::debug("- 单元数量: {}", exdyna.hexahedron_elements.size()); // 替换 std::cout
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
				spdlog::debug("为实例 {} 添加节点集 {} (包含 {} 个节点)", 
						instance_name, set_name, set_nodes.size()); // 替换 std::cout
			} 
			else {
				// 没有指定实例，视为装配级别的全局节点集 (旧行为保持兼容)
				std::map<std::string, std::vector<std::size_t>> empty_map;
				empty_map["__GLOBAL__"] = set_nodes; // 使用特殊键表示全局节点
				currentAssemblyPtr->node_sets[set_name] = empty_map;
				spdlog::debug("添加装配级全局节点集 {} (包含 {} 个节点)", 
						set_name, set_nodes.size()); // 替换 std::cout
			}
		}
		return true;
	}

	bool DataIn::BOUNDARY() {
		// 检查是否包含 op=NEW 参数
		bool op_new = (str.find("OP=NEW") != std::string::npos);
		bool is_velocity = (str.find("TYPE=VELOCITY") != std::string::npos);
		
		// 如果有op=NEW，设置对应类型的重置标志（不清除边界条件）
		if (op_new) {
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				if (is_velocity) {
					exdyna.steps.back().resetVelBoundary = true;
					spdlog::debug("步骤 {} 设置速度边界重置标志 (op=NEW)", exdyna.steps.back().name);
				} else {
					exdyna.steps.back().resetSpcBoundary = true;
					spdlog::debug("步骤 {} 设置位移边界重置标志 (op=NEW)", exdyna.steps.back().name);
				}
			} else {
				if (!exdyna.steps.empty()) {
					if (is_velocity) {
						exdyna.steps.front().resetVelBoundary = true;
						spdlog::debug("Initial步骤设置速度边界重置标志 (op=NEW)");
					} else {
						exdyna.steps.front().resetSpcBoundary = true;
						spdlog::debug("Initial步骤设置位移边界重置标志 (op=NEW)");
					}
				}
			}
		}
		
		// 检查下一行是否为另一个关键字
		std::streampos currentPos = fin.tellg(); // 保存当前位置
		std::string nextLine;
		std::getline(fin, nextLine);
		bool isEmptyBoundaryCommand = false;
		
		// 检查下一行是否为关键字或者结束
		if (nextLine.empty() || nextLine[0] == '*') {
			isEmptyBoundaryCommand = true;
			spdlog::debug("检测到空的边界条件命令: {}", str);
		}
		
		// 返回到当前命令行结束位置
		fin.seekg(currentPos);
		
		// 如果不是空命令，则读取并处理边界条件数据
		if (!isEmptyBoundaryCommand) {
			if (is_velocity) {
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
		// 检查是否包含 op=NEW 参数
		bool op_new = (str.find("OP=NEW") != std::string::npos);
		
		// 如果有op=NEW，设置对应类型的重置标志（不清除载荷数据）
		if (op_new) {
			// 如果在步骤定义中且有步骤数据，则设置重置标志
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 只设置重置标志，不清除数据
				exdyna.steps.back().resetDload = true;
				spdlog::debug("步骤 {} 设置分布载荷重置标志 (op=NEW)", exdyna.steps.back().name);
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
			spdlog::debug("位移边界条件已清除，无新添加的约束"); // 替换 std::cout
			return;
		}

		spdlog::debug("开始处理位移边界条件数据..."); // 替换 std::cout
		
		// 读取单点约束数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
            spdlog::debug("  解析位移边界行: {}", str); // 替换 std::cout
			
			// 创建临时约束数据 - 修改为使用std::string存储节点集名称
			std::pair<std::string, std::array<std::size_t, 3>> boundary_data;
			std::istringstream iss(str);
			std::string token;

			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				// 使用预定义的节点集合
				std::getline(iss, token, ',');
				spdlog::debug("  使用节点集: {}", token); // 替换 std::cout
				// 保存集合名称，不再查找集合中的节点
				boundary_data.first = token;
			}
			else {
				// 单个节点
				std::getline(iss, token, ',');
				spdlog::debug("  使用单个节点: {}", token); // 替换 std::cout
				// 使用#前缀表示单个节点
				boundary_data.first = "#" + token;
			}

			// 处理特殊约束类型
			if (str.find("PINNED") != std::string::npos || str.find("ENCASTRE") != std::string::npos) {
				spdlog::debug("  使用特殊约束类型: {}", (str.find("PINNED") != std::string::npos ? "PINNED" : "ENCASTRE")); // 替换 std::cout
				boundary_data.second.at(0) = 0;  // 起始自由度
				boundary_data.second.at(1) = 2;  // 结束自由度
				boundary_data.second.at(2) = 1;  // 增量
			}
			else if (str.find("XSYMM") != std::string::npos) {
				// X对称：约束X方向位移
				spdlog::debug("  使用特殊约束类型: XSYMM (X对称)"); // 替换 std::cout
				boundary_data.second.at(0) = 0;  // X方向 (第1个自由度)
				boundary_data.second.at(1) = 0;  // 仅X方向
				boundary_data.second.at(2) = 1;  // 增量
			}
			else if (str.find("YSYMM") != std::string::npos) {
				// Y对称：约束Y方向位移
				spdlog::debug("  使用特殊约束类型: YSYMM (Y对称)"); // 替换 std::cout
				boundary_data.second.at(0) = 1;  // Y方向 (第2个自由度)
				boundary_data.second.at(1) = 1;  // 仅Y方向
				boundary_data.second.at(2) = 1;  // 增量
			}
			else if (str.find("ZSYMM") != std::string::npos) {
				// Z对称：约束Z方向位移
				spdlog::debug("  使用特殊约束类型: ZSYMM (Z对称)"); // 替换 std::cout
				boundary_data.second.at(0) = 2;  // Z方向 (第3个自由度)
				boundary_data.second.at(1) = 2;  // 仅Z方向
				boundary_data.second.at(2) = 1;  // 增量
			}
			else {
				// 读取自由度约束定义
				std::getline(iss, token, ',');
				if (token.empty()) {
					// 如果没有指定自由度，默认约束所有自由度
					spdlog::debug("  未指定自由度，默认约束所有自由度 (1-3)"); // 替换 std::cout
					boundary_data.second.at(0) = 0;
					boundary_data.second.at(1) = 2;
					boundary_data.second.at(2) = 1;
				} else {
					spdlog::debug("  自由度起始值: {}", token); // 替换 std::cout
					boundary_data.second.at(0) = std::stoi(token) - 1;  // 起始自由度
					
					std::getline(iss, token, ',');
					if (token.empty()) {
						// 如果没有指定结束自由度，则设置为与起始自由度相同
						boundary_data.second.at(1) = boundary_data.second.at(0);
						spdlog::debug("  自由度结束值: 与起始值相同"); // 替换 std::cout
					}
					else {
						boundary_data.second.at(1) = std::stoi(token) - 1;  // 结束自由度
						spdlog::debug("  自由度结束值: {}", token); // 替换 std::cout
					}
					
					// 读取增量
					std::getline(iss, token, ',');
					if (token.empty()) {
						boundary_data.second.at(2) = 1;  // 默认增量为1
						spdlog::debug("  增量: 默认为1"); // 替换 std::cout
					} else {
						boundary_data.second.at(2) = std::stoi(token);
						spdlog::debug("  增量: {}", token); // 替换 std::cout
					}
				}
			}
            
			// 添加到当前步骤的约束列表中
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 添加到当前解析中的步骤
				exdyna.steps.back().boundary.spc_nodes.push_back(boundary_data);
				spdlog::debug("在步骤 {} 添加固定约束", exdyna.steps.back().name); // 替换 std::cout
			} else {
				// 如果不在Step定义中，则添加到Initial步骤，而不是全局约束
				if (!exdyna.steps.empty()) {
					exdyna.steps.front().boundary.spc_nodes.push_back(boundary_data);
					 spdlog::debug("在Initial步骤添加固定约束"); // 替换 std::cout
				}
			}
		}
		
		spdlog::debug("位移边界条件处理完成"); // 替换 std::cout
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
				spdlog::debug("引用实例 {} 中的", instance_name); // 替换 std::cout
			}
			
			// 处理节点集合或单个节点
			if (is_set) {
				// 使用预定义的节点集合
				std::string set_name = node_ref;
				// 构建完整的集合名称（包含实例前缀，如果有）
				if (!instance_name.empty()) {
					set_name = instance_name + "." + node_ref;
				}
				
				// 保存集合名称，而不是查找节点索引
				exdyna.ini_vel_generation.back().first = set_name;
				spdlog::debug("为节点集 {} 设置初始速度", set_name); // 替换 std::cout
			}
			else {
				// 单个节点，直接保存完整编号作为特殊的"集合名称"
				exdyna.ini_vel_generation.back().first = "#" + node_ref; // 使用#前缀标记它是单个节点ID
				spdlog::debug("为节点 {} 设置初始速度", node_ref); // 替换 std::cout
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
					spdlog::debug("  - 自由度: {}, 速度: {}", (exdyna.ini_vel_generation.back().second.first + 1), exdyna.ini_vel_generation.back().second.second); // 替换 std::cout
				} else {
					// 如果没有指定速度值，默认为0
					exdyna.ini_vel_generation.back().second.second = 0.0;
					spdlog::debug("  - 自由度: {}, 速度: 0.0 (默认)", (exdyna.ini_vel_generation.back().second.first + 1)); // 替换 std::cout
				}
			}
			catch (const std::exception& e) {
				spdlog::warn("解析速度数据时出错: {}", e.what()); // 替换 std::cerr
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
			spdlog::debug("速度边界条件已清除，无新添加的约束"); // 替换 std::cout
			return;
		}
		
		spdlog::debug("开始处理速度边界条件数据..."); // 替换 std::cout
		
		// 读取速度边界条件数据
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			spdlog::debug("  解析速度边界行: {}", str); // 替换 std::cout
			
			std::istringstream iss(str);
			std::string token;

			// 创建临时结构存储速度边界条件 - 修改为使用std::string存储节点集名称
			std::pair<std::string, std::pair<std::array<std::size_t, 2>, Types::Real>> boundary_data;
			boundary_data.second.first[0] = 0;
			boundary_data.second.first[1] = 0;
			boundary_data.second.second = 0.0;
			
			// 处理节点集合或单个节点
			if (str.find("SET-") != std::string::npos) {
				std::getline(iss, token, ',');
				spdlog::debug("  使用节点集: {}", token); // 替换 std::cout
				// 保存集合名称，不再查找集合中的节点
				boundary_data.first = token;
			}
			else {
				std::getline(iss, token, ',');
				 spdlog::debug("  使用单个节点: {}", token); // 替换 std::cout
				// 使用#前缀表示单个节点
				boundary_data.first = "#" + token;
			}

			// 读取自由度信息
			std::getline(iss, token, ',');
			if (token.empty()) {
				spdlog::warn("  警告: 缺少自由度信息"); // 替换 std::cerr
				continue; // 跳过这一行
			}
			
			spdlog::debug("  自由度起始值: {}", token); // 替换 std::cout
			boundary_data.second.first[0] = std::stoi(token) - 1;
			
			std::getline(iss, token, ',');
			if (token.empty()) {
				boundary_data.second.first[1] = boundary_data.second.first[0];
				spdlog::debug("  自由度结束值: 与起始值相同"); // 替换 std::cout
			}
			else {
				boundary_data.second.first[1] = std::stoi(token) - 1;
				spdlog::debug("  自由度结束值: {}", token); // 替换 std::cout
			}

			// 读取速度值
			if (std::getline(iss, token, ',')) {
				boundary_data.second.second = convertString<Types::Real>(token);
				 spdlog::debug("  速度值: {}", boundary_data.second.second); // 替换 std::cout
			} else {
				spdlog::debug("  默认速度值: 0.0"); // 替换 std::cout
			}
            
			// 添加到对应的数据结构中
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 添加到当前步骤的约束列表
				exdyna.steps.back().boundary.vel_nodes.push_back(boundary_data);
				spdlog::debug("在步骤 {} 添加速度约束", exdyna.steps.back().name); // 替换 std::cout
			} else {
				// 如果不在Step定义中，则添加到Initial步骤，而不是全局约束
				if (!exdyna.steps.empty()) {
					exdyna.steps.front().boundary.vel_nodes.push_back(boundary_data);
					spdlog::debug("在Initial步骤添加速度约束"); // 替换 std::cout
				}
			}
		}
		
		spdlog::debug("速度边界条件处理完成"); // 替换 std::cout
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
			spdlog::debug("为步骤 {} 设置重力加速度: {} 方向: ({}, {}, {})", // 替换 std::cout
				exdyna.steps.back().name, value, direction_x, direction_y, direction_z);
		} else {
			// 否则存储在全局变量中（不应该发生，但作为后备）
			spdlog::warn("在步骤外定义重力，这可能导致不一致的行为"); // 替换 std::cout
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
				spdlog::debug("为实例 {} 添加单元集 {} (包含 {} 个单元)", 
						instance_name, set_name, element_indices.size()); // 替换 std::cout
			} 
			else {
				// 没有指定实例，视为装配级别的全局单元集 (旧行为保持兼容)
				std::map<std::string, std::vector<std::size_t>> empty_map;
				empty_map["__GLOBAL__"] = element_indices; // 使用特殊键表示全局单元
				currentAssemblyPtr->element_sets[set_name] = empty_map;
				spdlog::debug("添加装配级全局单元集 {} (包含 {} 个单元)", 
						set_name, element_indices.size()); // 替换 std::cout
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
			spdlog::warn("无法解析表面名称"); // 替换 std::cerr
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
			spdlog::debug("在实例 {} 中定义表面: {}", instance_name, surf_name); // 替换 std::cout
		} else {
			spdlog::debug("定义表面: {}", surf_name); // 替换 std::cout
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
			
			spdlog::debug("  - 表面组成: 单元集={}, 面类型={}", full_ele_set_name, surf_face_type); // 替换 std::cout
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
		// 检查是否包含 op=NEW 参数
		bool op_new = (str.find("OP=NEW") != std::string::npos);
		
		// 如果有op=NEW，设置对应类型的重置标志（不清除载荷数据）
		if (op_new) {
			// 如果在步骤定义中且有步骤数据，则设置重置标志
			if (currentState == ParseState::STEP_DEFINITION && !exdyna.steps.empty()) {
				// 只设置重置标志，不清除数据
				exdyna.steps.back().resetDsload = true;
				spdlog::debug("步骤 {} 设置分布面载荷重置标志 (op=NEW)", exdyna.steps.back().name);
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
				spdlog::debug("为步骤 {} 添加分布面载荷: 表面={} 类型={} 值={}", 
					exdyna.steps.back().name, surf_name, load_type, load_value);
			} else {
				// 否则存储在全局变量中（不应该发生，但作为后备）
				spdlog::warn("在步骤外定义分布面载荷，这可能导致不一致的行为"); 
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
				spdlog::debug("从inp文件中读取到输出间隔: {}", exdyna.time_interval); // 替换 std::cout
				
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
					spdlog::debug("从inp文件的*Output设置中读取到number interval={}, 设置time_interval={}", number_interval, exdyna.time_interval); // 替换 std::cout
				}
			} catch (const std::exception& e) {
				spdlog::error("解析OUTPUT的number interval参数时出错: {}", e.what()); // 替换 std::cerr
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
		spdlog::debug("=== 开始转换数据到exDyna3D结构 ==="); // 替换 std::cout
		
		// 统计数量
		size_t totalNodes = 0;
		size_t totalElements = 0;
		
		// 1. 统计所有Part的节点和单元数量
		spdlog::debug("处理的部件数量: {}", inp_data.parts.size()); // 替换 std::cout
		for (const auto& part : inp_data.parts) {
			totalNodes += part.nodes.size();
			totalElements += part.elements.size();
			 spdlog::debug("- 部件 '{}': {} 个节点, {} 个单元", part.name, part.nodes.size(), part.elements.size()); // 替换 std::cout
		}
		
		// 如果没有Part数据（即inp_data未被填充），则直接返回，不进行转换
		if (inp_data.parts.empty()) {
			spdlog::warn("没有找到有效的Part数据，跳过转换。"); // 替换 std::cout
			return;
		}
		
		// 2. 为exDyna3D分配空间
		exdyna.vertices.clear();
		exdyna.hexahedron_elements.clear();
		exdyna.vertices.reserve(totalNodes);
		exdyna.hexahedron_elements.reserve(totalElements);
		
		// 3. 处理装配中的实例
		if (!inp_data.assemblies.empty()) {
			spdlog::debug("处理的装配数量: {}", inp_data.assemblies.size()); // 替换 std::cout
			
			for (auto& assembly : inp_data.assemblies) {
				spdlog::debug("- 处理装配 '{}' 及其实例", assembly.name); // 替换 std::cout
				
				size_t globalNodeIndex = 0;
				size_t globalElementIndex = 0;
				
				// 处理实例
				for (auto& instance : assembly.instances) {
					spdlog::debug("  - 实例 '{}' 引用部件 '{}'", instance.name, instance.part_name); // 替换 std::cout
					
					// 找到对应的部件
					auto partIt = std::find_if(inp_data.parts.begin(), inp_data.parts.end(), 
						[&](const InpPart& p) { return p.name == instance.part_name; });
					
					if (partIt != inp_data.parts.end()) {
						// 记录当前实例节点的起始索引
						instance.node_start_index = globalNodeIndex;
						instance.element_start_index = globalElementIndex;
						
						spdlog::debug("    节点起始索引: {}, 单元起始索引: {}", instance.node_start_index, instance.element_start_index); // 替换 std::cout
						
						// 检查是否需要应用平移变换
						bool hasTranslation = false;
						if (instance.translation[0] != 0.0 || instance.translation[1] != 0.0 || instance.translation[2] != 0.0) {
							hasTranslation = true;
							spdlog::debug("    应用平移变换: ({}, {}, {})", 
								instance.translation[0], instance.translation[1], instance.translation[2]);
						}
						
						// 添加节点到exdyna，应用平移变换
						for (const auto& node : partIt->nodes) {
							Types::Point<3> transformedNode = node;
							
							// 应用平移变换
							if (hasTranslation) {
								transformedNode[0] += instance.translation[0];
								transformedNode[1] += instance.translation[1];
								transformedNode[2] += instance.translation[2];
							}
							
							exdyna.vertices.push_back(transformedNode);
							globalNodeIndex++;
						}
						
						// 添加单元到exdyna
						for (const auto& element : partIt->elements) {
							exdyna.hexahedron_elements.emplace_back();
							auto& newElement = exdyna.hexahedron_elements.back();
							
							// 设置单元的物理ID和材料ID（默认为1）
							newElement.set_PID(1);
							newElement.set_MID(1);
							
							// 复制单元的节点索引，并调整为全局索引
							std::vector<Types::Vertex_index> globalNodes;
							globalNodes.reserve(element.size());
							
							for (const auto& localIdx : element) {
								globalNodes.push_back(localIdx + instance.node_start_index);
							}
							
							// 设置调整后的节点索引
							newElement.set_verticesIndex(globalNodes);
							globalElementIndex++;
						}
					} else {
						spdlog::warn("  警告: 找不到部件 '{}'", instance.part_name); // 替换 std::cerr
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
					spdlog::debug("  添加装配级节点集: {} 包含 {} 个节点", setName, globalNodeSet.size()); // 替换 std::cout
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
					spdlog::debug("  添加装配级单元集: {} 包含 {} 个单元", setName, globalElemSet.size()); // 替换 std::cout
				}
			}
		} else {
			// 如果没有装配数据，则直接添加所有Part
			spdlog::debug("未找到装配数据，将直接转换所有部件"); // 替换 std::cout
			
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
		}
		
		// 6. 处理材料属性
		if (!inp_data.materials.empty()) {
			spdlog::debug("处理材料数据: {} 个材料", inp_data.materials.size()); // 替换 std::cout
			
			// 这里只取第一个材料（假设单一材料）
			if (!inp_data.materials.empty()) {
				const auto& material = inp_data.materials[0];
				exdyna.mMatElastic.rho = material.density;
				exdyna.mMatElastic.E = material.E;
				exdyna.mMatElastic.v = material.v;
				spdlog::debug("设置材料属性: 密度={} E={} v={}", material.density, material.E, material.v); // 替换 std::cout
			}
		}
		
		spdlog::debug("数据转换完成: {} 个节点, {} 个单元", exdyna.vertices.size(), exdyna.hexahedron_elements.size()); // 替换 std::cout
		spdlog::debug("=== 转换完成 ===\n"); // 替换 std::cout
	}
}

