#pragma once
#include <string>
#include <fstream>
#include <map>
#include "exDyna3D.h"


namespace EnSC {
	class exDyna3D;
	class DataIn {
	public:
		DataIn(exDyna3D& p_exdyna);
		void read_inp(std::string fileName);
	private:
		exDyna3D& exdyna;
		bool NODE();
		bool ELEMENT();
		bool MATERIAL();
		bool SET_NODE_LIST();
		bool BOUNDARY();
		bool INITIAL_CONDITIONS();
		bool DLOAD();
		bool DSLOAD();
		bool TIME();
		bool AMP();
		bool SET_ELE_LIST();
		bool SURFACE();
		bool BULK_VISCOSITY();
		void SPC_NODE();
		void INITIAL_VELOCITY();
		void VELOCITY_NODE();
		void GRAV(std::string& amp_name);
		std::string str;
		std::ifstream fin;
		std::map<std::string, std::function<bool(DataIn*)>> k_func;
	};
}