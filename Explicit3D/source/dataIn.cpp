#include"../include/dataIn.h"
#include <algorithm>
#include <cctype>
#include <unordered_map>
#include <string>
#include <cctype>

namespace EnSC {
	template<typename T>
	inline T convertString(const std::string& str);

	template<>
	inline float convertString<float>(const std::string& str) {
		return std::stof(str);
	}

	template<>
	inline double convertString<double>(const std::string& str) {
		return std::stod(str);
	}

	inline static void toUpperCase(std::string& str) {
		std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
	}

	inline static void removeSpaces(std::string& str) {
		str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
	}

	DataIn::DataIn(exDyna3D& p_exdyna)
		:exdyna(p_exdyna) {
		k_func["*NODE"] = &DataIn::NODE;
		k_func["*ELEMENT"] = &DataIn::ELEMENT;
		k_func["*MATERIAL"] = &DataIn::MATERIAL;
		k_func["*BOUNDARY"] = &DataIn::BOUNDARY;
		k_func["*INITIAL CONDITION"] = &DataIn::INITIAL_CONDITIONS;
		k_func["*NSET"] = &DataIn::SET_NODE_LIST;
		k_func["*DLOAD"] = &DataIn::DLOAD;
		k_func["*DSLOAD"] = &DataIn::DSLOAD;
		k_func["*DYNAMIC"] = &DataIn::TIME;
		k_func["*AMPLITUDE"] = &DataIn::AMP;
		k_func["*ELSET"] = &DataIn::SET_ELE_LIST;
		k_func["*SURFACE"] = &DataIn::SURFACE;
		k_func["*BULK VISCOSITY"] = &DataIn::BULK_VISCOSITY;
	}

	void DataIn::read_inp(std::string fileName) {
		fin.open(fileName);
		if (!fin.is_open()) {
			std::cout << "failed to open the file!" << std::endl;
			exit(0);
		}
		else {
			std::cout << "Open " << fileName << " successfully!" << std::endl;
		}


		while (!fin.eof()) {
			if (fin.peek() == '*') {
				std::getline(fin, str);
				if (str.length() > 1) {
					char secondChar = str[1];
					if (secondChar != '*') {
						toUpperCase(str);
						for (const auto& pair : k_func) {
							const std::string& key = pair.first;
							// 检查 str 是否包含 key 作为子字符串
							if (str.find(key) != std::string::npos) {
								k_func[key](this);
								break;
							}
						}
					}
				}
			}
			else {
				std::getline(fin, str);
			}
		}

		fin.close();
		std::cout << "Read " << fileName << " Sucessfully!" << std::endl;
	}

	bool DataIn::NODE() {
		while (fin.peek() == '*') {
			std::getline(fin, str);
		}

		Types::Point<3> x;
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.empty()) {
				continue;
			}

			std::istringstream iss(str);
			std::string token;

			// Get nodeId
			std::getline(iss, token, ',');

			// Get x[0]
			std::getline(iss, token, ',');
			x[0] = convertString<Types::Real>(token);

			// Get x[1]
			std::getline(iss, token, ',');
			x[1] = convertString<Types::Real>(token);

			// Get x[2]
			std::getline(iss, token, ',');
			x[2] = convertString<Types::Real>(token);


			exdyna.vertices.emplace_back(x);
		}
		return true;
	}

	bool DataIn::MATERIAL() {
		Types::Real density = 0.0;
		Types::Real elasticModulus = 0.0;
		Types::Real poissonRatio = 0.0;
		std::string materialName;

		auto toUpperCase = [](std::string& str) {
			std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
			};

		while (std::getline(fin, str)) {
			toUpperCase(str);
			if (str.find("*MATERIAL") != std::string::npos) {
				auto pos = str.find("NAME=");
				if (pos != std::string::npos) {
					materialName = str.substr(pos + 5);
				}
			}
			else if (str.find("*DENSITY") != std::string::npos) {
				std::getline(fin, str);
				std::istringstream densityStream(str);
				densityStream >> density;
			}
			else if (str.find("*ELASTIC") != std::string::npos) {
				std::getline(fin, str);
				std::istringstream iss(str);
				std::string token;
				std::getline(iss, token, ',');
				elasticModulus = convertString<Types::Real>(token);
				std::getline(iss, token, ',');
				poissonRatio = convertString<Types::Real>(token);
			}
			else if (str.find("*") == 0 && str.find("*MATERIAL") == std::string::npos) {
				break;
			}
		}

		exdyna.mMatElastic.rho = density;
		exdyna.mMatElastic.E = elasticModulus;
		exdyna.mMatElastic.v = poissonRatio;
		exdyna.mMatElastic.update();

		return true;
	}

	bool DataIn::ELEMENT() {
		while (fin.peek() == '*') {
			std::getline(fin, str);
		}

		int PID = 1;
		int MID = 1;
		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);

			// 跳过空行
			if (str.empty()) continue;

			std::istringstream iss(str);
			std::string token;

			// 读取索引
			std::getline(iss, token, ',');

			// 读取节点索引
			std::vector<Types::Vertex_index> verticesIndex;
			while (std::getline(iss, token, ',')) {
				verticesIndex.push_back(std::stoul(token));
			}

			// 需要修改
			for (auto& idx : verticesIndex) {
				idx -= 1;
			}

			exdyna.hexahedron_elements.emplace_back();
			exdyna.hexahedron_elements.back().set_PID(PID);
			exdyna.hexahedron_elements.back().set_MID(MID);
			exdyna.hexahedron_elements.back().set_verticesIndex(verticesIndex);
		}
		return true;
	}

	bool DataIn::SET_NODE_LIST() {
		while (fin.peek() == '*')
			std::getline(fin, str);

		std::string::size_type pos = str.find("NSET=");
		std::string set_name;
		pos += 5; // Move past "NSET="
		while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
			set_name += str[pos++];
		}

		bool is_generate = (str.find("GENERATE") != std::string::npos);
		std::vector<std::size_t> set_nodes;

		std::getline(fin, str);

		std::istringstream iss(str);
		std::size_t value;
		while (iss >> value) {
			set_nodes.push_back(value);
			if (iss.peek() == ',') {
				iss.ignore();
			}
		}

		if (is_generate && set_nodes.size() == 3) {
			std::size_t start = set_nodes[0];
			std::size_t end = set_nodes[1];
			std::size_t step = set_nodes[2];
			set_nodes.clear();
			for (std::size_t i = start; i <= end; i += step) {
				set_nodes.push_back(i);
			}
		}

		for (auto& idx : set_nodes) {
			idx -= 1;
		}

		exdyna.map_set_node_list[set_name] = set_nodes;
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

	bool DataIn::INITIAL_CONDITIONS() {
		if (str.find("TYPE") != std::string::npos) {
			if (str.find("VELOCITY") != std::string::npos) {
				INITIAL_VELOCITY();
			}
		}
		return true;
	}

	bool DataIn::DLOAD() {
		std::string amp_name;
		std::string::size_type pos = str.find("AMPLITUDE=");
		if (pos != std::string::npos) {
			pos += 10;
			while (pos < str.size() && str[pos] != ',' && str[pos] != ' ') {
				amp_name += str[pos++];
			}
		}
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			if (str.find("GRAV") != std::string::npos) {
				GRAV(amp_name);
			}
		}

		return true;
	}

	bool DataIn::TIME() {
		while (fin.peek() == '*') {
			std::getline(fin, str);
		}

		while (fin.peek() != '*' && fin.peek() != EOF) {
			std::getline(fin, str);
			std::istringstream iss(str);
			std::string token;
			std::getline(iss, token, ',');
			std::getline(iss, token, ',');

			exdyna.totalTime = convertString<Types::Real>(token);
		}
		return true;
	}

	bool DataIn::BULK_VISCOSITY() {
		std::getline(fin, str);
		std::istringstream iss(str);
		std::string token;

		std::getline(iss, token, ',');
		exdyna.Cvisl = convertString<Types::Real>(token);
		std::getline(iss, token, ',');
		exdyna.Cvisq = convertString<Types::Real>(token);
		return true;
	}

	void DataIn::SPC_NODE() {
		while (fin.peek() == '*')
			std::getline(fin, str);

		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			exdyna.boundary_spc_node.emplace_back();
			std::istringstream iss(str);
			std::string token;

			//std::vector<std::pair<std::vector<std::size_t>, std::array<std::size_t, 3>>> boundary_spc_node;

			if (str.find("SET-") != std::string::npos) {
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().first = exdyna.map_set_node_list[token];

			}
			else {
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().first.push_back(std::stoi(token));
			}

			if (str.find("PINNED") != std::string::npos || str.find("ENCASTRE") != std::string::npos) {
				exdyna.boundary_spc_node.back().second.at(0) = 0;
				exdyna.boundary_spc_node.back().second.at(1) = 2;
				exdyna.boundary_spc_node.back().second.at(2) = 1;
			}
			else {
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().second.at(0) = std::stoi(token) - 1;
				std::getline(iss, token, ',');
				if (token.empty()) {
					exdyna.boundary_spc_node.back().second.at(1) = exdyna.boundary_spc_node.back().second.at(0);
				}
				else {
					exdyna.boundary_spc_node.back().second.at(1) = std::stoi(token) - 1;
				}
				std::getline(iss, token, ',');
				exdyna.boundary_spc_node.back().second.at(2) = std::stoi(token);

			}
		}
	}

	void DataIn::INITIAL_VELOCITY() {
		while (fin.peek() == '*')
			std::getline(fin, str);
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			exdyna.ini_vel_generation.emplace_back();
			std::istringstream iss(str);
			std::string token;


			if (str.find("SET-") != std::string::npos) {
				std::getline(iss, token, ',');
				exdyna.ini_vel_generation.back().first = exdyna.map_set_node_list[token];

			}
			else {
				std::getline(iss, token, ',');
				exdyna.ini_vel_generation.back().first.push_back(std::stoi(token) - 1);
			}

			std::getline(iss, token, ',');
			exdyna.ini_vel_generation.back().second.first = std::stoi(token) - 1;
			std::getline(iss, token, ',');
			exdyna.ini_vel_generation.back().second.second = convertString<Types::Real>(token);

		}
	}

	void DataIn::VELOCITY_NODE() {
		while (fin.peek() == '*')
			std::getline(fin, str);
		while (fin.peek() != '*' && !fin.eof()) {
			std::getline(fin, str);
			toUpperCase(str);
			std::istringstream iss(str);
			std::string token;

			std::pair<std::array<std::size_t, 2>, Types::Real>temp_pair;
			temp_pair.first[0] = 0;
			temp_pair.first[1] = 0;
			temp_pair.second = 0.0;
			if (str.find("SET-") != std::string::npos) {
				std::getline(iss, token, ',');
				exdyna.map_set_node_list[token];
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

