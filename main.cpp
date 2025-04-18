#include "exDyna3D.h"
#include "logger.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <cstdlib>
#include "spdlog/spdlog.h"

using namespace EnSC;
//主函数main

int main() {
	// 初始化日志系统
	Logger::init();

	spdlog::info("EnSCDynamics 显式动力学求解器启动");
	
	exDyna3D exdynarun;
	DataIn dataIn(exdynarun);
	//get k file name
	std::string str, filename;
	std::ifstream fin;
	fin.open("project.txt");
	if (!fin.is_open()) {
		spdlog::error("无法打开配置文件: project.txt");
		std::cout << "Can not open \"project.txt\"" << std::endl;
		exit(-1);
	}
	
	// 跳过注释行，读取inp文件名
	bool found_filename = false;
	while (std::getline(fin, str) && !found_filename) {
		// 忽略空行
		if (str.empty()) continue;
		// 忽略注释行
		if (str[0] == '#') continue;
		// 找到非注释行，非空行，视为文件名
		filename = str;
		found_filename = true;
	}
	
	fin.close();
	
	if (!found_filename) {
		spdlog::error("无法从project.txt中找到有效的inp文件名");
		std::cout << "无法从project.txt中找到有效的inp文件名" << std::endl;
		exit(-1);
	}
	
	filename = "inp/" + filename;
	spdlog::info("加载输入文件: {}", filename);
	dataIn.read_inp(filename);
	// Start the timer
	auto start = std::chrono::high_resolution_clock::now();

	// Run the exdynarun
	spdlog::info("开始运行求解器...");
	exdynarun.run();
	spdlog::info("求解器运行完成");

	// Stop the timer
	auto end = std::chrono::high_resolution_clock::now();

	// Calculate the duration
	std::chrono::duration<double> elapsed = end - start;
	spdlog::info("运行时间: {}秒", elapsed.count());
	std::cout << "exdynarun.run() took " << elapsed.count() << " seconds." << std::endl;
	
	// 关闭日志系统
	Logger::shutdown();

	system("pause");
	return 0;
}
