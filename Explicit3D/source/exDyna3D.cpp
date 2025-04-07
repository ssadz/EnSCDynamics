#include"../include/exDyna3D.h"
#include<iostream>
#include <unordered_set>
#include<omp.h>
#include <set>   // Include the header for std::set
#include <vector>
#include <utility> // For std::move
#include <numeric>

#include <random> // 用于生成随机数
#include <set>    // 用于存储不重复的随机块索引
#include <cmath>    // <<< 新增
#include <iomanip>  // <<< 新增 (用于设置输出精度)


namespace EnSC {
	using namespace Types;
	using namespace Eigen;
	void Material_elastic::update() {
		G = E / ((one + v) * two);
		K = E / (one - two * v) / three;
		lambda = two * G * v / (one - two * v);
		WOS = std::sqrt((lambda + two * G) / rho);
	}

	exDyna3D::exDyna3D() :
		nEle(0),
		gauss_intg_1d_1(1),
		feValuesHex(gauss_intg_1d_1, vertices),
		time((Real)0.0),
		totalTime((Real)0.0),
		dt_i((Real)0.0),
		dt_i_1((Real)0.0),
		factor_timeStep((Real)0.71),
		c_dr((Real)0.0125),
		c_cr((Real)0.03),
		Cvisl((Real)0.06),
		Cvisq((Real)1.2),
		time_interval(one),
		time_output(zero) {
		std::get<0>(gravity) = false;
	}

	exDyna3D::~exDyna3D() {}

	void exDyna3D::run() {
		init();
		get_fsiSph_virtualParticles_and_vel(2, 1.0/10.0);
		while (time < totalTime) {
			update_minVertex_perMaterial_and_dt();
			computeSate();
		}

	}

	// 在 exDyna3D.cpp 中
	void exDyna3D::init() {

		static bool first_flag = true;

		if (first_flag == true) {

			first_flag = false;

		}

		else {

			return;

		}

		this->read_project_txt();



		this->init_data();

		this->compute_invMassMatrix();

		this->apply_initial_condition();

		this->update_minVertex_perMaterial();

		this->apply_boundary_condition_vec();

		this->add_inForce_to_rhs();

		this->apply_external_node_force();

		this->apply_boundary_condition_a();

		this->compute_a();

		this->update_velocity();

		this->solution_a.setZero();

	}

	void exDyna3D::init_data() {
		static bool have_run = false;
		if (have_run) {
			return;
		}
		std::size_t nDofs = 3 * vertices.size();
		nEle = hexahedron_elements.size();
		inv_mass_matrix.resize(nDofs);
		solution_a.resize(nDofs);
		solution_rf.resize(nDofs);
		solution_v.resize(nDofs);
		solution_u.resize(nDofs);
		solution_du.resize(nDofs);
		solution_f.resize(nDofs);
		dtStable.resize(nEle);
		von_mises.resize(nEle);

		all_sigma.resize(nEle);
		all_Volume.resize(nEle);
		all_QIA.resize(nEle);
		all_vel_ele.resize(nEle);
		all_F.resize(nEle);
		all_Length.resize(nEle);
		all_VolRate.resize(nEle);
		all_xyz_matrix.resize(nEle);
		all_BT.resize(nEle);
		all_bT.resize(nEle);
		all_J.resize(nEle);

		inv_mass_matrix.setZero();
		solution_a.setZero();
		solution_v.setZero();
		solution_u.setZero();
		solution_rf.setZero();
		solution_du.setZero();
		solution_f.setZero();
		von_mises.setZero();
		dtStable.setZero();

		for (std::size_t i = 0; i < nEle; ++i) {
			std::fill(all_sigma[i].begin(), all_sigma[i].end(), 0.0);
			all_Volume[i] = 0.0;
			std::fill(all_QIA[i].begin(), all_QIA[i].end(), 0.0);
			all_vel_ele[i].setZero();
			// F 通常在 reinit 中设置为 I + du*BT，所以这里 Zero() 或 Identity() 可能都行，
			// 但 Zero() 更符合原始构造函数行为。检查 reinit_some_CellData 确认。
			// reinit_some_CellData 中 F = I + u_ele * BT，所以 Zero() 在这里是合适的
			all_F[i].setZero();
			all_Length[i] = 0.0;
			all_VolRate[i] = 0.0;
			all_xyz_matrix[i].setZero();
			all_BT[i].setZero();
			all_bT[i].setZero();
		}

		feValuesHex.init(EleType::hexN8);
		have_run = true;
	}

	void exDyna3D::output_results() {
		bool write_judger = (time >= time_output);
		bool end_flag = false;
		if (time >= totalTime) {
			end_flag = true;
		}
		if (write_judger || end_flag) {
			// It might be better to increment time_output only if a write actually happens,
			// but following the original logic here:
			time_output += time_interval;

			// Use a static counter for consistent file numbering across calls
			static int count = 0;

			// Create DataOut object
			DataOut dataOut(vertices, hexahedron_elements, time);

			// --- Prepare and Add FEM Data ---
			this->computeVonMises(); // Compute necessary derived data

			// Add node data
			dataOut.add_node_data(DataOutType::DataType::vector, &solution_a, "a");       // 加速度
			dataOut.add_node_data(DataOutType::DataType::vector, &solution_f, "force");   // 力 (Nodal forces calculated in FEM part)
			dataOut.add_node_data(DataOutType::DataType::vector, &solution_v, "v");       // 速度
			dataOut.add_node_data(DataOutType::DataType::vector, &solution_u, "u");       // 位移

			// Add mass matrix data (as nodal data)
			// Note: VTK typically visualizes scalar/vector/tensor per point/cell.
			// Representing a diagonal mass matrix might be better as a scalar field.
			Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> mass_matrix = inv_mass_matrix.cwiseInverse();
			// Consider adding just one component if it's lumped mass per node
			// dataOut.add_node_data(DataOutType::DataType::scalar, &mass_matrix, "mass");
			// Or if you really need 3 identical values per node:
			Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> mass_vector(mass_matrix.size() * 3);
			for (Eigen::Index i = 0; i < mass_matrix.size(); ++i) {
				mass_vector(3 * i + 0) = mass_matrix(i);
				mass_vector(3 * i + 1) = mass_matrix(i);
				mass_vector(3 * i + 2) = mass_matrix(i);
			}
			// If using the vector approach:
			// dataOut.add_node_data(DataOutType::DataType::vector, &mass_vector, "mass_vec");
			// For simplicity, let's assume scalar mass per node makes sense for visualization:
			Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> scalar_mass(vertices.size());
			for (size_t i = 0; i < vertices.size(); ++i) {
				// Assuming inv_mass_matrix stores x,y,z mass components consecutively
				// and they are the same for lumped mass. Take the first one.
				if (inv_mass_matrix(3 * i) != 0) // Avoid division by zero if mass is zero
					scalar_mass(i) = 1.0 / inv_mass_matrix(3 * i);
				else
					scalar_mass(i) = 0.0;
			}
			dataOut.add_node_data(DataOutType::DataType::scalar, &scalar_mass, "mass");


			// Add element data (Von Mises)
			dataOut.add_ele_data(DataOutType::DataType::scalar, &von_mises, "vonMises");

			// Prepare and Add Stress Tensor Data
			Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> stressVec(6 * nEle);
			for (std::size_t iEle = 0; iEle < nEle; ++iEle) {
				const auto& sigma = all_sigma[iEle];
				stressVec[6 * iEle + 0] = sigma[0]; // Sxx
				stressVec[6 * iEle + 1] = sigma[1]; // Syy
				stressVec[6 * iEle + 2] = sigma[2]; // Szz
				stressVec[6 * iEle + 3] = sigma[3]; // Sxy
				stressVec[6 * iEle + 4] = sigma[4]; // Sxz NOTE: Check order convention (was Syz in original?)
				stressVec[6 * iEle + 5] = sigma[5]; // Syz NOTE: Check order convention (was Szx in original?)
				// VTK symmetric tensor order is often xx, yy, zz, xy, yz, xz
				// Let's adjust to VTK common order if sigma stores it differently
				// Assuming sigma is [xx, yy, zz, xy, xz, yz]
				// If VTK wants [xx, yy, zz, xy, yz, xz]
				// stressVec[6 * iEle + 0] = sigma[0]; // xx
				// stressVec[6 * iEle + 1] = sigma[1]; // yy
				// stressVec[6 * iEle + 2] = sigma[2]; // zz
				// stressVec[6 * iEle + 3] = sigma[3]; // xy
				// stressVec[6 * iEle + 4] = sigma[5]; // yz <- Corrected index?
				// stressVec[6 * iEle + 5] = sigma[4]; // xz <- Corrected index?
				// **Action Required:** Verify the order of components stored in `allCellData_hexahedron[iEle].sigma`
				// and ensure it matches the order expected by VTK or your post-processor for symmetric tensors.
				// The code below assumes sigma stores [xx, yy, zz, xy, xz, yz] and VTK wants [xx, yy, zz, xy, yz, xz]
				stressVec[6 * iEle + 0] = sigma[0]; // xx
				stressVec[6 * iEle + 1] = sigma[1]; // yy
				stressVec[6 * iEle + 2] = sigma[2]; // zz
				stressVec[6 * iEle + 3] = sigma[3]; // xy
				stressVec[6 * iEle + 4] = sigma[5]; // yz
				stressVec[6 * iEle + 5] = sigma[4]; // xz
			}
			dataOut.add_ele_data(DataOutType::DataType::symmetric_tensor, &stressVec, "Stress");


			// Prepare and Add Strain Tensor Data
			Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> strainVec(6 * nEle);
			for (std::size_t iEle = 0; iEle < nEle; ++iEle) {
				const auto& F = all_F[iEle]; // 3x3 Deformation Gradient
				// Calculate Green-Lagrange strain E = 0.5 * (F^T * F - I)
				// Or engineering strain (small strain assumption) e = 0.5 * (grad(u) + grad(u)^T)
				// The original code calculates something like engineering strain but uses F instead of grad(u)
				// Let's calculate Green-Lagrange Strain which is more appropriate for F
				Eigen::Matrix<Types::Real, 3, 3> FTF = F.transpose() * F;
				Eigen::Matrix<Types::Real, 3, 3> E_mat = 0.5 * (FTF - Eigen::Matrix<Types::Real, 3, 3>::Identity());
				// Store in Voigt notation (xx, yy, zz, xy, yz, xz for VTK) - Note: Green-Lagrange has factor of 2 for shear
				strainVec[6 * iEle + 0] = E_mat(0, 0); // Exx
				strainVec[6 * iEle + 1] = E_mat(1, 1); // Eyy
				strainVec[6 * iEle + 2] = E_mat(2, 2); // Ezz
				strainVec[6 * iEle + 3] = E_mat(0, 1); // Exy (Note: This is E_xy, not gamma_xy=2*E_xy)
				strainVec[6 * iEle + 4] = E_mat(1, 2); // Eyz
				strainVec[6 * iEle + 5] = E_mat(0, 2); // Exz
			}
			dataOut.add_ele_data(DataOutType::DataType::symmetric_tensor, &strainVec, "Green-Lagrange Strain");


			// --- Write FEM Results ---
			const std::string femFileName = "solution_" + std::to_string(count) + ".vtu";
			dataOut.write_vtu(femFileName);
			// --- End FEM Output ---


			// --- Write Virtual Particle Results ---
			// Call the method added to DataOut, passing the current count and fsi_share_data
			dataOut.write_vtu_virtualParticles(count, this->fsi_share_data);
			// --- End Virtual Particle Output ---


			// Increment counter for the next output step
			count++;
		}
	}

	void exDyna3D::apply_initial_condition() {
		if (!ini_vel_generation.empty()) {
			for (auto iter = ini_vel_generation.begin(); iter != ini_vel_generation.end(); iter++) {
				for (auto& node_id : iter->first) {
					const auto dof0 = node_id * 3;
					const auto offset = iter->second.first;
					const auto value = iter->second.second;
					solution_v[dof0 + offset] = value;
				}
			}
		}
	}

	void exDyna3D::apply_boundary_condition_vec() {
		// Apply velocity boundary conditions
		if (!boundary_vel_node.empty()) {
			// Parallelize the outer loop over boundary conditions
#pragma omp parallel for
			for (std::size_t iter_idx = 0; iter_idx < boundary_vel_node.size(); ++iter_idx) {
				const auto& iter = boundary_vel_node[iter_idx];
				const auto& node_ids = iter.first;
				const auto& index_range = iter.second.first;
				const auto value = iter.second.second;

				// Loop over node IDs
				for (const auto& node_id : node_ids) {
					const auto dof0 = node_id * 3;
					// Loop over degrees of freedom
					for (auto i = index_range[0]; i <= index_range[1]; ++i) {
						solution_v[dof0 + i] = value;
					}
				}
			}
		}

		// Apply fixed boundary conditions
		if (!boundary_spc_node.empty()) {
			// Parallelize the outer loop over boundary conditions
#pragma omp parallel for
			for (std::size_t iter_idx = 0; iter_idx < boundary_spc_node.size(); ++iter_idx) {
				const auto& iter = boundary_spc_node[iter_idx];
				const auto& node_ids = iter.first;
				const auto& index_range = std::array<std::size_t, 2>{ iter.second[0], iter.second[1] };
				const auto condition = iter.second[2];

				// Loop over node IDs
				for (const auto& node_id : node_ids) {
					const auto dof0 = node_id * 3;
					// Loop over degrees of freedom
					for (auto i = index_range[0]; i <= index_range[1]; ++i) {
						if (condition == 1) {
							solution_v[dof0 + i] = 0.0;
						}
					}
				}
			}
		}
	}

	void exDyna3D::apply_boundary_condition_a() {
		solution_rf.setZero();

		// Apply velocity boundary conditions
		if (!boundary_vel_node.empty()) {
			// Parallelize the outer loop over boundary conditions
#pragma omp parallel for
			for (std::size_t iter_idx = 0; iter_idx < boundary_vel_node.size(); ++iter_idx) {
				const auto& iter = boundary_vel_node[iter_idx];
				const auto& node_ids = iter.first;
				const auto& index_range = iter.second.first;

				// Loop over node IDs
				for (const auto& node_id : node_ids) {
					const auto dof0 = node_id * 3;
					// Loop over degrees of freedom
					for (auto i = index_range[0]; i <= index_range[1]; ++i) {
						solution_a[dof0 + i] = 0.0;
					}
				}
			}
		}

		// Apply fixed boundary conditions
		if (!boundary_spc_node.empty()) {
			// Parallelize the outer loop over boundary conditions
#pragma omp parallel for
			for (std::size_t iter_idx = 0; iter_idx < boundary_spc_node.size(); ++iter_idx) {
				const auto& iter = boundary_spc_node[iter_idx];
				const auto& node_ids = iter.first;
				const auto index_start = iter.second[0];
				const auto index_end = iter.second[1];
				const auto condition = iter.second[2];

				// Loop over node IDs
				for (const auto& node_id : node_ids) {
					const auto dof0 = node_id * 3;
					// Loop over degrees of freedom
					for (auto i = index_start; i <= index_end; ++i) {
						if (condition == 1) {
							auto idx = dof0 + i;
							double sol_a = solution_a[idx];
							if (sol_a != 0.0) {
								solution_rf[idx] = -sol_a;
								solution_a[idx] = 0.0;
							}
						}
					}
				}
			}
		}
	}

	void exDyna3D::compute_invMassMatrix() {
		for (int iEle = 0; iEle < nEle; ++iEle) {
			using namespace Eigen;
			using Real = Real;
			feValuesHex.reinit(&hexahedron_elements[iEle]);
			all_Volume[iEle] = feValuesHex.get_JxW(0);

			auto& ele = hexahedron_elements[iEle];
			const auto& Ele_vertices = ele.get_verticesIndex();
			Real& rho = mMatElastic.rho;
			Real nodeMass = rho * all_Volume[iEle] / eight;

			all_BT[iEle] = feValuesHex.get_shape_derivatives_value_real_matrix(0).transpose(); // 新增

			for (int i = 0; i < 8; ++i) {
				auto& iVertex = Ele_vertices[i];
				const std::size_t dof = 3 * iVertex;
				inv_mass_matrix[dof] += nodeMass;
				inv_mass_matrix[dof + 1] += nodeMass;
				inv_mass_matrix[dof + 2] += nodeMass;
			}
		}
		inv_mass_matrix = inv_mass_matrix.cwiseInverse();
	}

	void exDyna3D::add_inForce_to_rhs() {
		// 定义总自由度数量
		const std::size_t nDofs = 3 * vertices.size();

		// 定义临时解向量，用于并行累加力
		// 为了避免竞态条件，使用多个临时向量
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a0;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a1;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a2;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a3;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a4;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a5;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a6;
		static Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a7;

		// 如果临时解向量大小不匹配，则调整大小
		if (nDofs != solution_a0.size()) {
			solution_a0.resize(nDofs);
			solution_a1.resize(nDofs);
			solution_a2.resize(nDofs);
			solution_a3.resize(nDofs);
			solution_a4.resize(nDofs);
			solution_a5.resize(nDofs);
			solution_a6.resize(nDofs);
			solution_a7.resize(nDofs);
		}

		// 将所有临时解向量初始化为零
		solution_a0.setZero();
		solution_a1.setZero();
		solution_a2.setZero();
		solution_a3.setZero();
		solution_a4.setZero();
		solution_a5.setZero();
		solution_a6.setZero();
		solution_a7.setZero();

		{
#pragma omp parallel for 
			for (std::size_t iEle = 0; iEle < nEle; ++iEle) {
				// 初始化应力张量为零
				Eigen::Matrix<Types::Real, 8, 3> f = Eigen::Matrix<Types::Real, 8, 3>::Zero();

				// 重新初始化单元数据
				reinit_some_CellData(iEle);

				// 计算应力响应
				// evaluateKirchhoffResponse(iEle, f); // 如果需要，可以启用
				evaluateJaumannResponse(iEle, f);

				// 计算沙漏力和体积粘度
				computeHourglassForce(iEle, f);
				computeVolumetricViscosity(iEle, f);

				// 计算稳定时间步长
				compute_dtStable(iEle);

				// 获取当前单元的节点索引
				auto& ele = hexahedron_elements[iEle];
				const auto& Ele_vertices = ele.get_verticesIndex();

				// 将应力张量的各分量映射为向量
				Eigen::Map<const Eigen::Vector<Types::Real, 8>> fx(&f(0, 0));
				Eigen::Map<const Eigen::Vector<Types::Real, 8>> fy(&f(0, 1));
				Eigen::Map<const Eigen::Vector<Types::Real, 8>> fz(&f(0, 2));

				// 遍历单元的所有8个节点
				for (std::size_t node = 0; node < 8; ++node) {
					std::size_t node_id = Ele_vertices[node];
					std::size_t dof0 = node_id * 3; // 三维情况下，每个节点有3个自由度

					// 根据节点索引，累加对应的力到临时解向量
					switch (node) {
					case 0:
						solution_a0[dof0] += fx[node];
						solution_a0[dof0 + 1] += fy[node];
						solution_a0[dof0 + 2] += fz[node];
						break;
					case 1:
						solution_a1[dof0] += fx[node];
						solution_a1[dof0 + 1] += fy[node];
						solution_a1[dof0 + 2] += fz[node];
						break;
					case 2:
						solution_a2[dof0] += fx[node];
						solution_a2[dof0 + 1] += fy[node];
						solution_a2[dof0 + 2] += fz[node];
						break;
					case 3:
						solution_a3[dof0] += fx[node];
						solution_a3[dof0 + 1] += fy[node];
						solution_a3[dof0 + 2] += fz[node];
						break;
					case 4:
						solution_a4[dof0] += fx[node];
						solution_a4[dof0 + 1] += fy[node];
						solution_a4[dof0 + 2] += fz[node];
						break;
					case 5:
						solution_a5[dof0] += fx[node];
						solution_a5[dof0 + 1] += fy[node];
						solution_a5[dof0 + 2] += fz[node];
						break;
					case 6:
						solution_a6[dof0] += fx[node];
						solution_a6[dof0 + 1] += fy[node];
						solution_a6[dof0 + 2] += fz[node];
						break;
					case 7:
						solution_a7[dof0] += fx[node];
						solution_a7[dof0 + 1] += fy[node];
						solution_a7[dof0 + 2] += fz[node];
						break;
					default:
						// 如果有更多节点，可以在此处添加处理
						break;
					}
				}
			}
		}
		// 将所有临时解向量累加到主解向量
		solution_a += (solution_a0 + solution_a1 + solution_a2 + solution_a3 +
			solution_a4 + solution_a5 + solution_a6 + solution_a7);
	}

	void exDyna3D::compute_a() {
		solution_a = solution_a.cwiseProduct(inv_mass_matrix);
	}

	void exDyna3D::update_velocity() {
		solution_v += (dt_i + dt_i_1) * solution_a / two;
		solution_du = dt_i * solution_v;
	}

	void exDyna3D::update_displacement() {
		solution_u += dt_i * solution_v;
		solution_a.setZero();
	}

	void exDyna3D::update_minVertex_perMaterial() {
#pragma omp parallel for
		for (int iEle = 0; iEle < nEle; ++iEle) {
			min_ele_edge(iEle);
		}
	}

	void exDyna3D::update_minVertex_perMaterial_and_dt() {
		update_minVertex_perMaterial();
		Real tmp_minVertexDist = std::numeric_limits<Real>::max();
#pragma omp parallel for reduction(min:tmp_minVertexDist)
		for (int iEle = 0; iEle < nEle; ++iEle) {
			const Real& local_minVertexDist = all_Length[iEle]; // 新增
			if (local_minVertexDist < tmp_minVertexDist) {
				tmp_minVertexDist = local_minVertexDist;
			}
		}
		//compute dt_i
		dt_i_1 = dt_i;
		const Real time_prev = time + dt_i_1;
		dt_i = tmp_minVertexDist / mMatElastic.WOS;
		Real tmp_dt = dtStable.minCoeff();
		dt_i = dt_i < tmp_dt ? dt_i : tmp_dt;
		dt_i *= factor_timeStep;
		if ((totalTime - time_prev) < dt_i && (totalTime - time_prev) > zero) {
			dt_i = totalTime - time_prev;
		}
	}

	void exDyna3D::move_mesh() {
		std::size_t nVertics = vertices.size();
#pragma omp parallel for
		for (int i = 0; i < nVertics; ++i) {
			vertices[i][0] += dt_i * solution_v[3 * i];
			vertices[i][1] += dt_i * solution_v[3 * i + 1];
			vertices[i][2] += dt_i * solution_v[3 * i + 2];
		}
	}

	void exDyna3D::read_project_txt() {
		std::string str, fileName = "project.txt";
		std::ifstream fin;
		fin.open(fileName);

		if (!fin.is_open()) {
			std::cout << "Can not open " << fileName << std::endl;
			exit(-1);
		}

		//k file name skip
		std::getline(fin, str);
		std::getline(fin, str);
		std::getline(fin, str);

		//c_dr
		std::getline(fin, str);
		fin >> c_dr;
		std::getline(fin, str);
		std::getline(fin, str);

		//c_cr
		std::getline(fin, str);
		fin >> c_cr;
		std::getline(fin, str);
		std::getline(fin, str);


		//factor of time step
		std::getline(fin, str);
		fin >> factor_timeStep;
		std::getline(fin, str);
		std::getline(fin, str);

		std::getline(fin, str);
		fin >> time_interval;
		if (time_interval > totalTime) {
			time_interval = totalTime;
		}

		fin.close();
	}

	void exDyna3D::computeSate() {
		time += dt_i_1;
		apply_boundary_condition_vec();
		add_inForce_to_rhs();
		apply_external_node_force();
		apply_boundary_condition_a();
		compute_a();
		output_results();
		printSate();
		update_velocity();
		update_displacement();
		move_mesh();

	}

	void exDyna3D::apply_external_node_force() {
		//gravity
		apply_gravity();
	}

	void exDyna3D::apply_gravity() {
		if (std::get<0>(gravity) == true) {
			Real amp_value = one;
			if (!std::get<1>(gravity).empty()) {
				const auto& amp_name = std::get<1>(gravity);
				amp_value = deal_amp(amp_name);
			}
			const auto& gravity_value = std::get<2>(gravity);
			const auto& direction_x = std::get<3>(gravity);
			const auto& direction_y = std::get<4>(gravity);
			const auto& direction_z = std::get<5>(gravity);

#pragma omp parallel for 
			for (int i = 0; i < solution_a.size(); i += 3) {
				solution_a[i] += amp_value * gravity_value * direction_x / inv_mass_matrix[i];
				solution_a[i + 1] += amp_value * gravity_value * direction_y / inv_mass_matrix[i + 1];
				solution_a[i + 2] += amp_value * gravity_value * direction_z / inv_mass_matrix[i + 2];
			}
		}
	}

	inline void exDyna3D::compute_dtStable(const std::size_t& iEle) {
		const Real& Len = all_Length[iEle]; // 新增
		const Real& VolRate = all_VolRate[iEle]; // 新增
		const Real& cd = mMatElastic.WOS;
		const Real xi = Cvisl - Cvisq * Cvisq * Len / cd * std::min(zero, VolRate);
		dtStable[iEle] = (std::sqrt(one + xi * xi) - xi) * Len / cd;
	}

	inline void exDyna3D::reinit_some_CellData(const std::size_t& iEle) {
		auto& ele = hexahedron_elements[iEle];
		const auto& Ele_vertices = ele.get_verticesIndex();

		Eigen::Matrix<Real, 3, 8> u_ele;
		u_ele.setZero();

		for (std::size_t i = 0; i < 8; ++i) {
			const auto node_id = Ele_vertices[i];
			const std::size_t dof0 = 3 * node_id;

			// 设置速度
			all_vel_ele[iEle](0, i) = solution_v[dof0];
			all_vel_ele[iEle](1, i) = solution_v[dof0 + 1];
			all_vel_ele[iEle](2, i) = solution_v[dof0 + 2];

			// 设置位移
			u_ele(0, i) = solution_u[dof0];
			u_ele(1, i) = solution_u[dof0 + 1];
			u_ele(2, i) = solution_u[dof0 + 2];

			// 设置 xyz_matrix
			const auto& tmpPoint = vertices[node_id];
			all_xyz_matrix[iEle](i, 0) = tmpPoint[0];
			all_xyz_matrix[iEle](i, 1) = tmpPoint[1];
			all_xyz_matrix[iEle](i, 2) = tmpPoint[2];
		}

		// 计算变形梯度 F = I + u_ele * BT
		all_F[iEle] = Eigen::Matrix<Real, 3, 3>::Identity() + u_ele * all_BT[iEle];

		all_J[iEle] = all_F[iEle].determinant();

		// 计算 bT = BT * F.inverse()
		all_bT[iEle] = all_BT[iEle] * all_F[iEle].inverse();

		// 映射 bT 的各列为向量
		Eigen::Map<const Eigen::Vector<Real, 8>> b0(&all_bT[iEle](0, 0));
		Eigen::Map<const Eigen::Vector<Real, 8>> b1(&all_bT[iEle](0, 1));
		Eigen::Map<const Eigen::Vector<Real, 8>> b2(&all_bT[iEle](0, 2));

		// 计算体积变化率 VolRate = vel_ele.row(0).dot(b0) + vel_ele.row(1).dot(b1) + vel_ele.row(2).dot(b2)
		all_VolRate[iEle] = all_vel_ele[iEle].row(0).dot(b0) + all_vel_ele[iEle].row(1).dot(b1) + all_vel_ele[iEle].row(2).dot(b2);
	}

	inline void exDyna3D::evaluateJaumannResponse(const std::size_t& iEle, Eigen::Matrix<Real, 8, 3>& f) {
		// 获取当前单元的数据
		Eigen::Map<const Eigen::Vector<Real, 8>> xI(&all_xyz_matrix[iEle](0, 0));
		Eigen::Map<const Eigen::Vector<Real, 8>> yI(&all_xyz_matrix[iEle](0, 1));
		Eigen::Map<const Eigen::Vector<Real, 8>> zI(&all_xyz_matrix[iEle](0, 2));


		// 计算速度梯度 L = vel_ele * bT (3x8 * 8x3 = 3x3)
		Eigen::Matrix<Real, 3, 3> L = all_vel_ele[iEle] * all_bT[iEle];

		// 计算应变增量 epsilonInc = [εxx, εyy, εzz, εxy, εxz, εyz]
		Eigen::Matrix<Real, 6, 1> epsilonInc;
		epsilonInc(0) = L(0, 0); // εxx
		epsilonInc(1) = L(1, 1); // εyy
		epsilonInc(2) = L(2, 2); // εzz
		epsilonInc(3) = (L(0, 1) + L(1, 0)) / two; // εxy
		epsilonInc(4) = (L(0, 2) + L(2, 0)) / two; // εxz
		epsilonInc(5) = (L(1, 2) + L(2, 1)) / two; // εyz

		epsilonInc *= dt_i_1; // 乘以时间步长

		// 计算旋转速度 omega = [ω01, ω02, ω12]
		Eigen::Matrix<Real, 3, 1> omega;
		omega(0) = (L(0, 1) - L(1, 0)) / two; // ω01
		omega(1) = (L(0, 2) - L(2, 0)) / two; // ω02
		omega(2) = (L(1, 2) - L(2, 1)) / two; // ω12

		// 获取材料参数
		const Real lambda = mMatElastic.lambda;
		const Real G = mMatElastic.G;

		// 计算Jaumann应力增量
		Eigen::Matrix<Real, 6, 1> Jaumann;
		Jaumann(0) = lambda * (epsilonInc(0) + epsilonInc(1) + epsilonInc(2)) + 2 * G * epsilonInc(0); // σxx
		Jaumann(1) = lambda * (epsilonInc(0) + epsilonInc(1) + epsilonInc(2)) + 2 * G * epsilonInc(1); // σyy
		Jaumann(2) = lambda * (epsilonInc(0) + epsilonInc(1) + epsilonInc(2)) + 2 * G * epsilonInc(2); // σzz
		Jaumann(3) = 2 * G * epsilonInc(3); // σxy
		Jaumann(4) = 2 * G * epsilonInc(4); // σxz
		Jaumann(5) = 2 * G * epsilonInc(5); // σyz

		// 映射sigma向量
		Eigen::Map<Eigen::Vector<Real, 6>> _sigma(&all_sigma[iEle][0]);

		// 计算应力增量 stressInc = [Δσxx, Δσyy, Δσzz, Δσxy, Δσxz, Δσyz]
		Eigen::Matrix<Real, 6, 1> stressInc;
		stressInc(0) = Jaumann(0) + 2 * _sigma(3) * omega(0) * dt_i_1 + 2 * _sigma(4) * omega(1) * dt_i_1; // Δσxx
		stressInc(1) = Jaumann(1) - 2 * _sigma(3) * omega(0) * dt_i_1 + 2 * _sigma(5) * omega(2) * dt_i_1; // Δσyy
		stressInc(2) = Jaumann(2) - 2 * _sigma(4) * omega(1) * dt_i_1 - 2 * _sigma(5) * omega(2) * dt_i_1; // Δσzz
		stressInc(3) = Jaumann(3) + (_sigma(1) - _sigma(0)) * omega(0) * dt_i_1 + _sigma(5) * omega(1) * dt_i_1; // Δσxy
		stressInc(4) = Jaumann(4) + (_sigma(2) - _sigma(0)) * omega(1) * dt_i_1 + _sigma(5) * omega(0) * dt_i_1; // Δσxz
		stressInc(5) = Jaumann(5) + (_sigma(2) - _sigma(1)) * omega(2) * dt_i_1 + _sigma(3) * omega(1) * dt_i_1 - _sigma(4) * omega(0) * dt_i_1; // Δσyz

		// 更新sigma
		_sigma += stressInc;
		const Real& Volume = all_Volume[iEle]; // 初始体积

		const Real& J = all_J[iEle];
		Real currentVolume = J * Volume; // 当前体积

		// 映射 bT 的各列为向量
		Eigen::Map<const Eigen::Vector<Real, 8>> b0(&all_bT[iEle](0, 0)); // 新增
		Eigen::Map<const Eigen::Vector<Real, 8>> b1(&all_bT[iEle](0, 1)); // 新增
		Eigen::Map<const Eigen::Vector<Real, 8>> b2(&all_bT[iEle](0, 2)); // 新增

		// 计算力的贡献 (使用当前体积 currentVolume)
		Eigen::Vector<Real, 8> fx_m = (b0 * _sigma(0) + b1 * _sigma(3) + b2 * _sigma(4)) * currentVolume; // 新增
		Eigen::Vector<Real, 8> fy_m = (b0 * _sigma(3) + b1 * _sigma(1) + b2 * _sigma(5)) * currentVolume; // 新增
		Eigen::Vector<Real, 8> fz_m = (b0 * _sigma(4) + b1 * _sigma(5) + b2 * _sigma(2)) * currentVolume; // 新增

		// 映射力矩阵 f
		Eigen::Map<Eigen::Vector<Real, 8>> fx(&f(0, 0)); // x方向力
		Eigen::Map<Eigen::Vector<Real, 8>> fy(&f(0, 1)); // y方向力
		Eigen::Map<Eigen::Vector<Real, 8>> fz(&f(0, 2)); // z方向力

		// 更新力
		fx -= fx_m;
		fy -= fy_m;
		fz -= fz_m;
	}

	inline void exDyna3D::computeHourglassForce(const std::size_t& iEle, Eigen::Matrix<Types::Real, 8, 3>& f) {
		// 沙漏控制常量
		const Types::Real C2 = 0.1;
		const Eigen::Matrix<Types::Real, 4, 8> TAU{
			{  1., -1.,  1., -1.,  1., -1.,  1., -1. },
			{ -1., -1.,  1.,  1.,  1.,  1., -1., -1. },
			{ -1.,  1.,  1., -1.,  1., -1., -1.,  1. },
			{ -1.,  1., -1.,  1.,  1., -1.,  1., -1. }
		};

		Types::Real Real_4divide3 = 4.0 / 3.0;
		Types::Real Real_1divideS8 = 1.0 / std::sqrt(8.0);

		// 定义变量
		Eigen::Matrix<Types::Real, 3, 4> HIA = Eigen::Matrix<Types::Real, 3, 4>::Zero();
		Eigen::Matrix<Types::Real, 8, 3> HGS = Eigen::Matrix<Types::Real, 8, 3>::Zero();
		Eigen::Matrix<Types::Real, 8, 3> HGD = Eigen::Matrix<Types::Real, 8, 3>::Zero();
		Eigen::Matrix<Types::Real, 4, 8> ATAU = Eigen::Matrix<Types::Real, 4, 8>::Zero();
		//Types::Real NX2, DR, CR, V;

		ATAU = TAU - TAU * all_xyz_matrix[iEle] * all_BT[iEle].transpose();
		HIA = all_vel_ele[iEle] * ATAU.transpose(); // 新增

		// 计算 NX2
		const auto NX2 = all_BT[iEle].squaredNorm(); // 新增


		// 获取单元体积和初始体积
		const Real& J = all_J[iEle];
		const Real& V0 = all_Volume[iEle];
		const Real V = J * V0;

		// 获取材料参数
		const Real& G = mMatElastic.G;
		const Real& K = mMatElastic.K;
		const Real& rho = mMatElastic.rho;
		const Real& WOS = mMatElastic.WOS;

		// 计算 DR 和 CR
		const auto DR = c_dr * (Real_4divide3 * G + K) * C2 * NX2 * V * Real_1divideS8 * dt_i_1;
		const auto CR = c_cr * rho * V0 / V * std::cbrt(V0 * V0) * WOS / four; // 使用 std::cbrt 替换 pow(V, 2.0/3.0) 可能略好


		// 更新 QIA (重要修改)
		for (int i = 0; i < 3; ++i) {
			for (int k = 0; k < 4; ++k) {
				int ik = 4 * i + k;
				// cellData.QIA[ik] += DR * HIA(i, k); // 删除
				all_QIA[iEle][ik] += DR * HIA(i, k); // 新增
			}
		}

		// 计算 HGS (重要修改)
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 8; ++j)
				for (int k = 0; k < 4; ++k)
				{
					int ik = 4 * i + k;
					// HGS(j, i) += cellData.QIA[ik] * ATAU(k, j); // 删除
					HGS(j, i) += all_QIA[iEle][ik] * ATAU(k, j); // 新增
				}

		// 计算 HGD (不变, 使用 HIA)
		HGD = ATAU.transpose() * HIA.transpose();
		HGD *= CR;

		// 合并 HGS 和 HGD (不变)
		HGS += HGD;

		// 将沙漏力添加到总力 f (不变)
		f -= HGS;
	}

	inline void exDyna3D::computeVolumetricViscosity(const std::size_t& iEle, Eigen::Matrix<Types::Real, 8, 3>& f) {

		const Real& V0 = all_Volume[iEle]; // 新增
		const Real& J = all_J[iEle];
		const Real V = J * V0; // 当前体积

		const Real rho_int = mMatElastic.rho; // 初始密度
		const Real c = mMatElastic.WOS; // 波速
		const Real Len = std::cbrt(V0);
		const Real& VolRate = all_VolRate[iEle]; // 新增
		const Real rho = rho_int / J; // 当前密度
		const Real p = rho * Len * (Cvisq * Cvisq * Len * VolRate * std::min(zero, VolRate) - Cvisl * c * VolRate);

		Eigen::Map<const Eigen::Vector<Real, 8>> b0(&all_bT[iEle](0, 0));
		Eigen::Map<const Eigen::Vector<Real, 8>> b1(&all_bT[iEle](0, 1));
		Eigen::Map<const Eigen::Vector<Real, 8>> b2(&all_bT[iEle](0, 2));

		// 力是 p * b * dV_current = p * b * V
		const Eigen::Vector<Real, 8> fx_v = p * b0 * V;
		const Eigen::Vector<Real, 8> fy_v = p * b1 * V;
		const Eigen::Vector<Real, 8> fz_v = p * b2 * V;

		Eigen::Map<Eigen::Vector<Real, 8>> fx(&f(0, 0));
		Eigen::Map<Eigen::Vector<Real, 8>> fy(&f(0, 1));
		Eigen::Map<Eigen::Vector<Real, 8>> fz(&f(0, 2));

		fx += fx_v;
		fy += fy_v;
		fz += fz_v;
	}

	void exDyna3D::computeVonMises() {
#pragma omp parallel for
		for (int iEle = 0; iEle < nEle; iEle++) {
			const auto& sigma = all_sigma[iEle];

			// 提取应力分量
			Real sigma_xx = sigma[0];
			Real sigma_yy = sigma[1];
			Real sigma_zz = sigma[2];
			Real sigma_xy = sigma[3];
			Real sigma_yz = sigma[4];
			Real sigma_zx = sigma[5];

			// 计算Von Mises应力
			Real term1 = (sigma_xx - sigma_yy) * (sigma_xx - sigma_yy);
			Real term2 = (sigma_yy - sigma_zz) * (sigma_yy - sigma_zz);
			Real term3 = (sigma_zz - sigma_xx) * (sigma_zz - sigma_xx);
			Real term4 = six * (sigma_xy * sigma_xy + sigma_yz * sigma_yz + sigma_zx * sigma_zx);

			Real vonMises = std::sqrt((term1 + term2 + term3 + term4) / two);
			von_mises[iEle] = vonMises;
		}
	}

	void exDyna3D::printSate() {
		static int result_Real = 0;
		if (result_Real % 5000 == 0) {
			printf("Soild Time: %.5e ; Time step: %d ; dt: %.5e\n", time, result_Real, dt_i);
		}
		result_Real++;
	}

	Types::Real exDyna3D::min_ele_edge(const std::size_t& iEle) {

		auto& ele = hexahedron_elements[iEle];
		const auto& Ele_vertices = ele.get_verticesIndex();
		Types::Real dx, dy, dz, minEleEdge, tmp;

		// 初始化为一个较大的值
		minEleEdge = std::numeric_limits<Types::Real>::max();

		// 计算底部面的边：0-1, 1-2, 2-3, 3-0
		for (int i = 0; i < 4; ++i)
		{
			int next_i = (i < 3) ? i + 1 : 0; // 确保循环到起点
			dx = vertices[Ele_vertices[next_i]][0] - vertices[Ele_vertices[i]][0];
			dy = vertices[Ele_vertices[next_i]][1] - vertices[Ele_vertices[i]][1];
			dz = vertices[Ele_vertices[next_i]][2] - vertices[Ele_vertices[i]][2];
			tmp = dx * dx + dy * dy + dz * dz;
			if (tmp < minEleEdge)
				minEleEdge = tmp;
		}

		// 计算顶部面的边：4-5, 5-6, 6-7, 7-4
		for (int i = 4; i < 8; ++i)
		{
			int next_i = (i < 7) ? i + 1 : 4; // 确保循环到起点
			dx = vertices[Ele_vertices[next_i]][0] - vertices[Ele_vertices[i]][0];
			dy = vertices[Ele_vertices[next_i]][1] - vertices[Ele_vertices[i]][1];
			dz = vertices[Ele_vertices[next_i]][2] - vertices[Ele_vertices[i]][2];
			tmp = dx * dx + dy * dy + dz * dz;
			if (tmp < minEleEdge)
				minEleEdge = tmp;
		}

		// 计算垂直边：0-4, 1-5, 2-6, 3-7
		for (int i = 0; i < 4; ++i)
		{
			dx = vertices[Ele_vertices[i + 4]][0] - vertices[Ele_vertices[i]][0];
			dy = vertices[Ele_vertices[i + 4]][1] - vertices[Ele_vertices[i]][1];
			dz = vertices[Ele_vertices[i + 4]][2] - vertices[Ele_vertices[i]][2];
			tmp = dx * dx + dy * dy + dz * dz;
			if (tmp < minEleEdge)
				minEleEdge = tmp;
		}

		minEleEdge = std::sqrt(minEleEdge);
		all_Length[iEle] = minEleEdge;
		return minEleEdge;
	}

	Real exDyna3D::deal_amp(const std::string& amp_name) {
		const auto& amp_vector = map_amp_list[amp_name];
		Real dt = dt_i;
		Real y = one;
		for (int i = 0; i < amp_vector.size() - 1; i++) {
			const Real x1 = amp_vector[i][0];
			const Real y1 = amp_vector[i][1];

			const Real x2 = amp_vector[i + 1][0];
			const Real y2 = amp_vector[i + 1][1];
			if (time > x2) {
				break;
			}
			else {
				const Real m = (y2 - y1) / (x2 - x1);
				const Real c = y1 - m * x1;
				Real x = time + dt / two;
				if ((x + dt / two) < x2) {
					y = m * x + c;
					break;
				}
				else {
					const Real DT1 = x2 - time;
					Real X = time + DT1 / two;
					Real Y = m * X + c;
					const Real A1 = Y * DT1;
					Real A2 = zero;
					for (int j = 2; j + i < amp_vector.size(); j++) {
						const Real xj = amp_vector[i + j][0];
						const Real yj = amp_vector[i + j][1];
						const Real xj_1 = amp_vector[i + j - 1][0];
						const Real yj_1 = amp_vector[i + j - 1][1];
						if ((x + dt / two) > xj) {

							A2 += (xj - xj_1) * (yj + yj_1) / two;
						}
						else {
							const Real M = (yj - yj_1) / (xj - xj_1);
							const Real C = yj_1 - m * xj_1;
							Real Y1 = M * (x + dt / two) + C;
							A2 += (x + dt / two - xj_1) * (yj_1 + Y1) / two;
							break;
						}
					}
					const Real x_end = (amp_vector.end() - 1)->at(0);
					if (x + dt / two > x_end) {
						//const Real y_end = (amp_vector.end() - 1)->at(1);
						A2 += (x + dt / two - x_end) * one;
					}

					y = (A1 + A2) / dt;
					break;
				}
			}
		}
		return y;
	}


    //生成虚粒子
	void exDyna3D::get_fsi_sph_Ele_face_index() {
		// 创建节点到包含该节点的单元数量的映射
		// Using std::map instead of std::unordered_map
		std::map<int, int> node_element_count;

		// 首先统计每个节点被多少个单元包含
		for (int iEle = 0; iEle < nEle; ++iEle)
		{
			const auto& ele = hexahedron_elements[iEle];
			const auto& eleVertices = ele.get_verticesIndex(); // Assuming this returns a container like std::vector<int>

			// 遍历单元的每个节点，增加计数
			for (int i = 0; i < 8; ++i)
			{
				int node_index = eleVertices[i];
				node_element_count[node_index]++; // std::map supports operator[] for insertion/access
			}
		}

		// 然后找出所有被少于8个单元包含的节点，这些就是边界节点
		// Using std::set instead of std::unordered_set
		std::set<int> usFsiSphNodes;
		// Range-based for loop is C++11
		for (const auto& node_count_pair : node_element_count) // Iterate over std::map pairs
		{
			// For std::map, the element is a std::pair<const Key, T>
			// node_count_pair.first is the key (node_index)
			// node_count_pair.second is the value (count)

			// 注意：对于六面体单元，内部节点通常被8个单元共享。
			// 边界节点被少于8个单元共享。角点可能被更少单元共享。
			// 这里假设与流体接触的面上的节点会被少于8个实体单元共享。
			// 这个判断条件可能需要根据具体情况调整（例如，如果实体模型本身就有内部空腔）。
			// *** 强烈建议后续改为使用用户定义的边界集 ***
			if (node_count_pair.second < 8)
			{
				usFsiSphNodes.insert(node_count_pair.first); // Insert the node index (key)
			}
		}


		// 识别与流体接触的单元面
		fsiBoundaryElementFaces.clear(); // 清空旧数据 (Assuming fsiBoundaryElementFaces is e.g. std::vector<std::pair<Types::Vertex_index, std::vector<int>>>)

		for (Types::Vertex_index iEle = 0; iEle < nEle; ++iEle)
		{
			const auto& ele = hexahedron_elements[iEle];
			const auto& eleVertices = ele.get_verticesIndex();
			std::vector<int> currentElementFaces; // 临时存储当前单元的边界`面

			// 遍历单元的6个面
			for (int iFace = 0; iFace < 6; ++iFace)
			{
				int fsiSphFaceJudger = 0; // 记录当前面上有多少个节点是边界节点
				const auto& faceIndices = ele.Ele_face_indices[iFace]; // 获取当前面的局部节点索引

				// 检查面上的4个节点是否都在边界节点集合中
				for (int i = 0; i < 4; ++i)
				{
					int node_index = eleVertices[faceIndices[i]]; // 获取全局节点索引
					// std::set also has the count method
					if (usFsiSphNodes.count(node_index))
					{
						++fsiSphFaceJudger;
					}
				}

				// 如果面的所有4个节点都在FSI边界上
				if (fsiSphFaceJudger == 4)
				{
					currentElementFaces.push_back(iFace); // 将边界`面索引加入临时列表
				}
			}
			// 如果当前单元找到了任何边界`面，则将其加入主列表
			if (!currentElementFaces.empty()) {
				// std::emplace_back and std::move are C++11 features
				fsiBoundaryElementFaces.emplace_back(iEle, std::move(currentElementFaces));
			}
		}
	}

	void exDyna3D::get_fsiSph_virtualParticles_and_vel(int nLayers, Types::Real virtualParticlesDist) {
		get_fsi_sph_Ele_face_index(); // 识别边界单元和面

		if (fsiBoundaryElementFaces.empty()) {
			std::cout << "警告：未找到 FSI 边界单元/面。" << std::endl;
			fsi_share_data.FSI_virtualParticles_coordinates.clear();
			fsi_share_data.FSI_virtualParticles_velocity.clear();
			fsi_share_data.FSI_virtualParticles_nodeForce.clear();
			map_eleIndex_virParticle_unit_data.clear();
			map_eleIndex_virParticles_index.clear();
			return;
		}

		// 调用修改后的 processInteractionElements，传递 nLayers 和 virtualParticlesDist
		processInteractionElements(nLayers, virtualParticlesDist);

		// resizeFsiSharedData 仍然需要，以调整力的向量大小
		resizeFsiSharedData();

		std::cout << "生成了 " << fsi_share_data.FSI_virtualParticles_coordinates.size()
				  << " 个最终虚拟粒子（沿法线生成，可能经过合并）。" << std::endl;
	}

	void exDyna3D::populateElementMatrices(int eleIdx, Eigen::Matrix<Types::Real, 8, 3>& xyzMatrix,
		Eigen::Matrix<Types::Real, 8, 3>& velMatrix) const {
		// 检查 eleIdx 的有效性
		if (eleIdx < 0 || eleIdx >= hexahedron_elements.size()) {
			throw std::out_of_range("Element index out of range in populateElementMatrices.");
		}

		const auto& ele = hexahedron_elements[eleIdx];
		const auto& verticesIdx = ele.get_verticesIndex(); // 获取单元的8个全局节点索引

		if (verticesIdx.size() != 8) {
			// 添加错误处理或断言
			throw std::runtime_error("Hexahedron element does not have 8 vertices.");
		}

		for (int i = 0; i < 8; ++i) { // 遍历单元的8个节点
			Types::Vertex_index globalNodeIdx = verticesIdx[i]; // 获取全局节点索引

			// 检查全局节点索引是否有效
			if (globalNodeIdx >= vertices.size()) {
				throw std::out_of_range("Global node index out of range for vertices vector.");
			}
			// 使用 .size() 检查比 >= 更安全
			if (3 * globalNodeIdx + 2 >= solution_v.size()) {
				throw std::out_of_range("Global node index out of range for solution_v vector.");
			}


			int dof = 3 * globalNodeIdx; // 计算该节点在全局自由度向量中的起始索引

			// 填充速度矩阵
			velMatrix(i, 0) = solution_v[dof];     // x方向速度
			velMatrix(i, 1) = solution_v[dof + 1]; // y方向速度
			velMatrix(i, 2) = solution_v[dof + 2]; // z方向速度

			// 填充坐标矩阵
			const auto& point = vertices[globalNodeIdx]; // 获取节点坐标
			xyzMatrix(i, 0) = point[0]; // x坐标
			xyzMatrix(i, 1) = point[1]; // y坐标
			xyzMatrix(i, 2) = point[2]; // z坐标
		}
	}

    static Eigen::Matrix<Types::Real, 3, 8> compute_shape_derivatives_at_point(const std::array<Types::Real, 3>& unitPoint) {
		Eigen::Matrix<Types::Real, 3, 8> derivatives_matrix_unit;
	    Element_HexN8 temp_ele; // 需要一个单元实例来调用方法
	    for (int iNode = 0; iNode < 8; ++iNode) {
         // 调用 Element_HexN8 中计算三个方向导数的函数
         std::array<Types::Real, 3> derivs = temp_ele.get_shapeFunctionDerivatives(iNode, unitPoint);
         derivatives_matrix_unit(0, iNode) = derivs[0]; // dNi/dxi
         derivatives_matrix_unit(1, iNode) = derivs[1]; // dNi/deta
         derivatives_matrix_unit(2, iNode) = derivs[2]; // dNi/dzeta
		 }
		 return derivatives_matrix_unit;
	}
    
	void exDyna3D::generateVirtualParticlesForFace(
		int faceIdx,                       // 当前处理的面索引 (0-5)
		int nLayers,                       // 要生成的粒子层数
		Types::Real virtualParticlesDist,  // 层之间的目标物理距离
		const Eigen::Matrix<Types::Real, 8, 3>& xyzMatrix, // 单元节点的当前物理坐标
		const Eigen::Matrix<Types::Real, 8, 3>& velMatrix, // 单元节点的当前速度
		std::vector<std::array<Types::Real, 3>>& tempElementCoords, // 输出: 生成的粒子坐标
		std::vector<std::array<Types::Real, 3>>& tempElementVels,   // 输出: 生成的粒子速度
		std::vector<std::array<Types::Real, 3>>& tempElementUnitCoords // 输出: 生成的粒子对应的单位坐标
	) {
		// --- 步骤 1: 确定面中心点和单位法线方向 ---
		std::array<Types::Real, 3> unitPointSurface; // 面中心的单位坐标
		Eigen::Matrix<Types::Real, 3, 1> n_unit_vec; // 单位坐标系中的法线向量 (指向单元内部)

		// 根据面索引 faceIdx 确定面中心坐标和法线方向
		switch (faceIdx)
		{
			case 0: unitPointSurface = {-1.0, 0.0, 0.0}; n_unit_vec << 1.0, 0.0, 0.0; break; // xi=-1 面, 法线 +xi
			case 1: unitPointSurface = { 1.0, 0.0, 0.0}; n_unit_vec << -1.0, 0.0, 0.0; break; // xi= 1 面, 法线 -xi
			case 2: unitPointSurface = { 0.0, -1.0, 0.0}; n_unit_vec << 0.0, 1.0, 0.0; break; // eta=-1 面, 法线 +eta
			case 3: unitPointSurface = { 0.0, 1.0, 0.0}; n_unit_vec << 0.0, -1.0, 0.0; break; // eta= 1 面, 法线 -eta
			case 4: unitPointSurface = { 0.0, 0.0, -1.0}; n_unit_vec << 0.0, 0.0, 1.0; break; // zeta=-1 面, 法线 +zeta
			case 5: unitPointSurface = { 0.0, 0.0, 1.0}; n_unit_vec << 0.0, 0.0, -1.0; break; // zeta= 1 面, 法线 -zeta
			default:
				// 增加健壮性：记录错误或抛出异常
				std::cerr << "错误: 在 generateVirtualParticlesForFace 中遇到无效的面索引 (" << faceIdx << ")" << std::endl;
				return; // 停止处理这个无效的面
		}

		// --- 步骤 2: 计算面中心点的雅可比矩阵 ---
		// 调用上面定义的辅助函数计算形函数导数
		Eigen::Matrix<Types::Real, 3, 8> shape_derivatives_unit_at_surface = compute_shape_derivatives_at_point(unitPointSurface);
		// 计算雅可比矩阵 J(3x3) = d(x,y,z)/d(xi,eta,zeta) = (dN/dXi)(3x8) * X(8x3)
		Eigen::Matrix<Types::Real, 3, 3> J = shape_derivatives_unit_at_surface * xyzMatrix;

		// --- 步骤 3: 计算尺度因子 s ---
		// 将单位法线映射到物理空间的法线方向向量
		Eigen::Matrix<Types::Real, 3, 1> n_real = J * n_unit_vec;
		// 计算物理空间法向量的模长的平方
		Types::Real s_squared = n_real.squaredNorm();
		Types::Real s = 0.0; // 初始化尺度因子
		// 定义一个小的容差值，避免除零或处理退化单元时出问题
		const Types::Real tolerance_squared = 1e-24; // 比较平方值以避免开方

		if (s_squared > tolerance_squared) {
			s = std::sqrt(s_squared); // 计算尺度因子 s = ||J * n_unit||
		} else {
			// 处理退化映射情况：无法可靠地确定步长
			// 记录警告信息，并跳过这个面的粒子生成
			std::cerr << "警告: 面 " << faceIdx << " 检测到退化映射或零尺度因子 (s^2 = " << s_squared
					<< ")。跳过此面的粒子生成。" << std::endl;
			return; // 退出此函数，不为这个面生成粒子
		}

		// --- 步骤 4: 计算单位坐标步长 ---
		// 根据期望的物理距离和尺度因子，计算在单位坐标系中需要的步长
		Types::Real unit_step = virtualParticlesDist / s;

		// --- 步骤 5: 生成粒子层 ---
		Element_HexN8 ele_for_shape_func; // 创建一个单元实例用于计算形函数值
		Eigen::Matrix<Types::Real, 1, 8> shapeFuncValue; // 用于存储形函数值的矩阵 N(1x8)

		// 从面上的点 (k=0) 开始，向内生成 nLayers 层粒子
		for (int k = 0; k < nLayers; ++k) {
			// 计算当前第 k 层的单位坐标点
			std::array<Types::Real, 3> unitPointLayer;
			unitPointLayer[0] = unitPointSurface[0] + k * unit_step * n_unit_vec(0);
			unitPointLayer[1] = unitPointSurface[1] + k * unit_step * n_unit_vec(1);
			unitPointLayer[2] = unitPointSurface[2] + k * unit_step * n_unit_vec(2);

			// 计算 unitPointLayer 处的形函数值 N
			for (int n = 0; n < 8; ++n) {
				// 使用单元实例获取节点 n 在 unitPointLayer 处的形函数值
				shapeFuncValue(0, n) = ele_for_shape_func.get_shapeFunctionValue(n, unitPointLayer);
			}

			// 插值计算物理坐标: realCoord(1x3) = N(1x8) * X(8x3)
			auto realCoordMat = shapeFuncValue * xyzMatrix;
			std::array<Types::Real, 3> coord = {realCoordMat(0, 0), realCoordMat(0, 1), realCoordMat(0, 2)};

			// 插值计算速度: velocity(1x3) = N(1x8) * V(8x3)
			auto velocityMat = shapeFuncValue * velMatrix;
			std::array<Types::Real, 3> vel = {velocityMat(0, 0), velocityMat(0, 1), velocityMat(0, 2)};

			// 存储计算得到的坐标、速度和对应的单位坐标
			tempElementCoords.push_back(coord);
			tempElementVels.push_back(vel);
			tempElementUnitCoords.push_back(unitPointLayer); // 存储用于生成的单位坐标
		}
	}

	void exDyna3D::setUnitPointForFace(int faceIdx, Types::Real u, Types::Real v, Types::Real offset,
		std::array<Types::Real, 3>& unitPoint) const {
		// 根据面的索引 (0-5)，将二维的单位坐标 (u, v) 和法向偏移 offset
		// 映射到三维的单位坐标 unitPoint = {xi, eta, zeta}
		// offset 是从面指向单元内部的距离（在单位坐标系下）
		switch (faceIdx)
		{
		case 0: unitPoint = { -1.0 + offset, u, v }; break;  // xi = -1 面 (对应节点 3,0,4,7) -> 法向是 +xi
		case 1: unitPoint = { 1.0 - offset, u, v }; break;   // xi = 1  面 (对应节点 1,2,6,5) -> 法向是 -xi
		case 2: unitPoint = { u, -1.0 + offset, v }; break;  // eta = -1 面 (对应节点 0,1,5,4) -> 法向是 +eta
		case 3: unitPoint = { u, 1.0 - offset, v }; break;   // eta = 1  面 (对应节点 2,3,7,6) -> 法向是 -eta
		case 4: unitPoint = { u, v, -1.0 + offset }; break;  // zeta = -1 面 (对应节点 1,0,3,2) -> 法向是 +zeta
		case 5: unitPoint = { u, v, 1.0 - offset }; break;   // zeta = 1  面 (对应节点 4,5,6,7) -> 法向是 -zeta
		default:
			// 无效的面索引，可以抛出异常或记录错误
			throw std::invalid_argument("Invalid face index (" + std::to_string(faceIdx) + ") provided to setUnitPointForFace.");
			break;
		}
	}

	void exDyna3D::resizeFsiSharedData() {
		// 获取当前生成的虚粒子数量
		int nParticles = fsi_share_data.FSI_virtualParticles_coordinates.size();

		// 调整存储虚粒子受力的向量大小
		// 通常在SPH计算完力之后，将力传递回来时才需要填充这个向量
		// 这里可以预先分配大小，并可能初始化为0
		try {
			// 只有在大小不匹配时才resize，可能更高效
			if (fsi_share_data.FSI_virtualParticles_nodeForce.size() != nParticles) {
				fsi_share_data.FSI_virtualParticles_nodeForce.resize(nParticles);
				// 可以选择将力初始化为零 (如果需要确保初始状态)
				// std::fill(fsi_share_data.FSI_virtualParticles_nodeForce.begin(),
				//           fsi_share_data.FSI_virtualParticles_nodeForce.end(),
				//           std::array<Types::Real, 3>{0.0, 0.0, 0.0});
			}
		}
		catch (const std::bad_alloc& e) {
			std::cerr << "Memory allocation failed when resizing FSI_virtualParticles_nodeForce: " << e.what() << std::endl;
			// 根据需要处理内存分配失败的情况
			throw; // 重新抛出异常，让上层知道失败了
		}
		catch (...) {
			std::cerr << "An unexpected error occurred during resizing FSI_virtualParticles_nodeForce." << std::endl;
			throw; // 重新抛出未知异常
		}
	}

	// ... existing code ...
    void exDyna3D::processInteractionElements(int nLayers, Types::Real virtualParticlesDist) {
        // --- 准备最终的全局列表 ---
        std::vector<std::array<Types::Real, 3>> finalGlobalCoords;
        std::vector<std::array<Types::Real, 3>> finalGlobalVels;

        // --- 处理单元前清除内部映射表 ---
        map_eleIndex_virParticle_unit_data.clear();
        map_eleIndex_virParticles_index.clear();

        // --- 遍历边界单元 ---
        for (const auto& elementFacePair : fsiBoundaryElementFaces) {
            int eleIdx = elementFacePair.first;
            const auto& faces = elementFacePair.second;
            int nFaces = faces.size();
            if (nFaces <= 0) continue;

            // --- 当前单元生成粒子的临时存储 ---
            std::vector<std::array<Types::Real, 3>> tempCoords;
            std::vector<std::array<Types::Real, 3>> tempVels;
            std::vector<std::array<Types::Real, 3>> tempUnitCoords; // 存储原始单位坐标

            // --- 为该单元生成所有潜在的粒子 ---
            Eigen::Matrix<Types::Real, 8, 3> xyzMatrix, velMatrix;
            populateElementMatrices(eleIdx, xyzMatrix, velMatrix);

            // 对该单元的每个边界`面，调用粒子生成函数
            for (int faceIdx : faces) {
                generateVirtualParticlesForFace(
                    faceIdx, 
                    nLayers,
                    virtualParticlesDist,  // 添加虚拟粒子距离参数
                    xyzMatrix, 
                    velMatrix,
                    tempCoords, 
                    tempVels, 
                    tempUnitCoords
                );
            }

            // --- 将 *所有* 生成的粒子添加到全局列表并更新映射表 ---
            int numParticles = tempCoords.size(); // 获取实际生成的粒子数
            size_t currentGlobalIndexStart = finalGlobalCoords.size();
            auto& elementUnitData = map_eleIndex_virParticle_unit_data[eleIdx];
            auto& elementGlobalIndices = map_eleIndex_virParticles_index[eleIdx];
            elementUnitData.clear();
            elementGlobalIndices.clear();
            
            // 预分配空间以提高效率
            elementUnitData.reserve(numParticles);
            elementGlobalIndices.reserve(numParticles);

            for (int k = 0; k < numParticles; ++k) {
                // 因为没有合并，所有粒子都是活动的，直接添加
                size_t globalIndex = currentGlobalIndexStart + k; // 计算最终的全局索引
                finalGlobalCoords.push_back(tempCoords[k]);
                finalGlobalVels.push_back(tempVels[k]);
                elementUnitData.push_back(tempUnitCoords[k]);
                elementGlobalIndices.push_back(globalIndex); // 存储对应的全局索引
            }
        }

        // --- 使用最终列表更新主要的 fsi_share_data ---
        fsi_share_data.FSI_virtualParticles_coordinates = std::move(finalGlobalCoords);
        fsi_share_data.FSI_virtualParticles_velocity = std::move(finalGlobalVels);
    } // <--- Adjusted this brace
// ... existing code ...



}