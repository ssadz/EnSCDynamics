#include "exDyna3D.h"
using namespace EnSC;
//主函数main

int main() {
	exDyna3D exdynarun;
	DataIn dataIn(exdynarun);
	//get k file name
	std::string str, filename;
	std::ifstream fin;
	fin.open("project.txt");
	if (!fin.is_open()) {
		std::cout << "Can not open " << "project.txt" << std::endl;
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
		std::cout << "无法从project.txt中找到有效的inp文件名" << std::endl;
		exit(-1);
	}
	
	filename = "inp/" + filename;
	dataIn.read_inp(filename);
	// Start the timer
	auto start = std::chrono::high_resolution_clock::now();

	// Run the exdynarun
	exdynarun.run();

	// Stop the timer
	auto end = std::chrono::high_resolution_clock::now();

	// Calculate the duration
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "exdynarun.run() took " << elapsed.count() << " seconds." << std::endl;
	system("pause");
	return 0;
}
