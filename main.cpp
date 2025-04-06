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
	std::getline(fin, str);
	fin >> filename;
	fin.close();
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
