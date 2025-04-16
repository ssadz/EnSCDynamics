#pragma once

// --- C++ 标准库头文件 ---
#include <vector>       // 包含 std::vector
#include <array>        // 包含 std::array
#include <map>          // 包含 std::map
#include <tuple>        // 包含 std::tuple
#include <string>       // 包含 std::string (可能被间接引入，但这里明确包含aa它)
#include <chrono>       // 包含时间相关 (如果需要计时功能)
#include <numeric>      // 包含 std::iota (如果需要 RCM 优化)
#include <algorithm>    // 包含 std::sort 等 (如果需要 RCM 优化)
#include <set>          // 包含 std::set (如果需要 RCM 或 FSI 边界处理)
#include <utility>      // 包含 std::pair, std::move

// --- Eigen 头文件 (通过 types.h 被引入) ---
// 通常不需要在这里显式包含 <Eigen/...>，而是直接使用的 Eigen 类型

// --- 项目内部头文件 ---
#include "types.h"          // 项目自定义类型 (包含 Eigen, Point 等)
#include "ElementHexN8.h"   // 六面体单元类
#include "gauss_integral_1d.h" // 高斯积分
#include "fevalues.h"       // 有限元值类
#include "dataOut.h"        // 数据输出
#include "dataIn.h"         // 数据输入
#include "FSI_share_data.h" // FSI 共享数据结构

namespace EnSC {

	// --- 前向声明 ---
	class DataIn;

	// --- 材料属性结构体 ---
	struct Material_elastic {
		void update(); // 更新派生材料属性的方法

		// 基本属性
		Types::Real E = static_cast<Types::Real>(2.06e11); // 杨氏模量
		Types::Real rho = static_cast<Types::Real>(7850.0); // 密度
		Types::Real v = static_cast<Types::Real>(0.3);   // 泊松比

		// 派生属性
		Types::Real G;      // 剪切模量 = E / (2 * (1 + v))
		Types::Real K;      // 体积模量 = E / (3 * (1 - 2 * v))
		Types::Real lambda; // 拉梅第一参数 = (E * v) / ((1 + v) * (1 - 2 * v))
		Types::Real WOS;    // 声速 (材料参数)
	};

	// --- 步骤数据结构 ---
	struct StepData {
		std::string name;                 // 步骤名称
		Types::Real timePeriod;           // 时间周期
		std::vector<std::pair<std::vector<std::size_t>, std::array<std::size_t, 3>>> boundary_spc_node; // 位移约束节点集合
		std::vector<std::pair<std::vector<std::size_t>, std::pair<std::array<std::size_t, 2>, Types::Real>>> boundary_vel_node; // 速度约束节点集合
		// 其他特定于步骤的数据可以在这里添加
	};

	// --- 单元数据结构体 (模板前向声明) ---
	// template<std::size_t n_integral_points> // 这个结构体定义似乎被移除 (改为 SoA)
	// struct CellDataHex;                     // 这个定义可能并不需要前向声明

	// --- 外部节点力结构体 (似乎未使用，可能需清理) ---
	// struct ExternalNodeForce {
	// 	std::vector<Types::Vertex_index> external_node_force_nodeID; // 节点ID
	// 	std::vector<int> external_node_force_direction;         // 作用方向
	// 	std::vector<Types::Real> external_node_force_value;         // 力的大小
	// };

	// --- 主求解器类 ---
	class exDyna3D {
	public:
		// --- 构造/析构函数 ---
		exDyna3D();
		~exDyna3D(); // 析构函数（可能需要释放资源，打印统计）

		// --- 友元类 (允许 DataIn 访问 protected/private 成员) ---
		friend class DataIn;

		// --- 主要公共接口 ---
		virtual void run();  // 执行模拟
		void init();         // 初始化模型
		
		// --- 时间步骤控制 ---
		void setCurrentStep(std::size_t stepIndex); // 设置当前时间步骤
		std::size_t getCurrentStep() const { return currentStepIndex; } // 获取当前时间步骤索引

		// --- FSI 共享数据 (允许外部访问/修改) ---
		FSI_share_data fsi_share_data;

		void calculate_time_interval(); // 计算输出时间间隔

	protected:
		// --- 初始化相关方法 ---
		virtual void init_data();           // 初始化数据结构，分配大小
		void read_project_txt();            // 读取配置文件 (project.txt)
		void compute_invMassMatrix();       // 计算质量矩阵逆 (对角质量)
		void apply_initial_condition();     // 施加初始条件 (如初始速度)

		// --- 主计算循环和状态更新 ---
		virtual void computeSate();         // 执行一个时间步的计算
		void add_inForce_to_rhs();          // 将内部力累加到右端项 (内部力)
		virtual void compute_a();           // 计算加速度 a = M^{-1} * F
		void update_velocity();             // 更新节点速度 (中心差分)
		void update_displacement();         // 更新节点位移 (中心差分)
		void move_mesh();                   // 基于位移/速度移动网格节点坐标

		// --- 单元相关计算 (在 add_inForce_to_rhs 调用) ---
		inline void reinit_some_CellData(const std::size_t& iEle); // 更新单元当前状态 (F, J, bT, VolRate 等)
		inline void evaluateJaumannResponse(const std::size_t& iEle, Eigen::Matrix<Types::Real, 8, 3>& f); // 计算 Jaumann 应力率更新及内部力数组
		inline void computeHourglassForce(const std::size_t& iEle, Eigen::Matrix<Types::Real, 8, 3>& f); // 计算沙漏力
		inline void computeVolumetricViscosity(const std::size_t& iEle, Eigen::Matrix<Types::Real, 8, 3>& f); // 计算体积粘性力

		// --- Helper functions for reinit_some_CellData ---
		void update_element_kinematics(const std::size_t& iEle);
		void compute_deformation_gradient(const std::size_t& iEle, const Eigen::Matrix<Types::Real, 3, 8>& u_ele); // Pass u_ele
		void compute_jacobian(const std::size_t& iEle);
		void compute_spatial_gradient(const std::size_t& iEle);
		void compute_volume_rate(const std::size_t& iEle);
		// --- End Helper functions ---

		// --- 边界条件相关方法 ---
		void apply_boundary_condition_vec(); // 施加 Dirichlet 边界条件 (速度/位移约束) 到速度数组
		void apply_boundary_condition_a();   // 施加 Dirichlet 边界条件到加速度数组 (约束残差力)
		void apply_external_node_force();    // 施加外力 (目前只有节点力)
		void apply_gravity();                // 施加重力

		// --- 时间步相关方法 ---
		void update_minVertex_perMaterial();        // 计算每种材料单元的最小边 (调用 min_ele_edge)
		void update_minVertex_perMaterial_and_dt(); // 计算并更新稳定时间步长 dt_i
		inline void compute_dtStable(const std::size_t& iEle); // 计算单个单元的粘性/声速稳定时间步
		Types::Real min_ele_edge(const std::size_t& iEle);     // 计算单元最小边长

		// --- FSI 相关方法 ---
		void get_fsi_sph_Ele_face_index();          // 获取 FSI 边界单元面索引
		void get_fsiSph_virtualParticles_and_vel(int nLayers, Types::Real virtualParticlesDist); // 生成 FSI 虚拟粒子
		void processInteractionElements(int nLayers, Types::Real virtualParticlesDist); // 处理交互单元及粒子（沿法线生成，含合并）
		void populateElementMatrices(int eleIdx, Eigen::Matrix<Types::Real, 8, 3>& xyzMatrix, Eigen::Matrix<Types::Real, 8, 3>& velMatrix) const; // 获取单元节点坐标和速度
		void apply_fsiSph_nodeForce();//把虚粒子的力应用到单元上
		void update_virParticles_coor_vel();//更新虚粒子坐标和速度


		// 新声明示例 (假设它仍在类中声明)
		void setUnitPointForFace(int faceIdx, Types::Real x, Types::Real y, Types::Real offset, std::array<Types::Real, 3>& unitPoint) const; // 设置虚拟粒子单位坐标
		void resizeFsiSharedData();                 // 调整 FSI 共享数据大小
		void generateVirtualParticlesForFace(
		int faceIdx,
		int nLayers,
		Types::Real virtualParticlesDist,
		const Eigen::Matrix<Types::Real, 8, 3>& xyzMatrix,
		const Eigen::Matrix<Types::Real, 8, 3>& velMatrix,
		std::vector<std::array<Types::Real, 3>>& tempElementCoords,
		std::vector<std::array<Types::Real, 3>>& tempElementVels,
		std::vector<std::array<Types::Real, 3>>& tempElementUnitCoords
	);

		// --- 输出和工具方法 ---
		void output_results();              // 输出结果文件
		void printSate();                   // 打印当前状态到终端
		void computeVonMises();             // 计算 Von Mises 应力
		Types::Real deal_amp(const std::string& amp_name); // 处理幅值曲线

		// --- 成员变量 ---
	protected: // 也可以是 private:

		// --- 几何和单元相关数据 ---
		Types::VerticesAll<3> vertices;         // 节点坐标数组 (共用)
		std::vector<Element_HexN8> hexahedron_elements; // 单元数组 (共用, 包含所有单元节点坐标)
		Types::Vertex_index nEle;               // 单元数量

		// --- 单元状态数据 (SoA, 按单元存储, 共用) ---
		std::vector<std::array<Types::Real, 6>> all_sigma; // 应力 [xx,yy,zz,xy,xz,yz] ? (确保顺序)
		std::vector<Types::Real> all_Volume;           // 单元参考体积 V0
		std::vector<std::array<Types::Real, 12>> all_QIA; // 沙漏应变率更新及内部力数组
		std::vector<Eigen::Matrix<Types::Real, 3, 8>> all_vel_ele; // 单元节点速度 (瞬时)
		std::vector<Eigen::Matrix<Types::Real, 3, 3>> all_F; // 单元力矩阵 F
		std::vector<Types::Real> all_Length;           // 单元长度 (瞬时)
		std::vector<Types::Real> all_VolRate;          // 体积变化 tr(L)
		std::vector<Eigen::Matrix<Types::Real, 8, 3>> all_xyz_matrix; // 单元节点当前坐标 (瞬时)
		std::vector<Eigen::Matrix<Types::Real, 8, 3>> all_BT; // 全单元力矩阵转置 (参考) B0^T
		std::vector<Eigen::Matrix<Types::Real, 8, 3>> all_bT; // 全单元力矩阵转置 (前一步) B^T
		std::vector<Types::Real> all_J;                // 单元力矩阵行列式 det(F) (瞬时)

		// --- 全解向量 (共用节点自由度存储, 共用) ---
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_u;  // 位移向量
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_v;  // 速度向量
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_a;  // 加速度向量 (也可以瞬时存储)
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_du; // 位移增量 (瞬时存储)
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_f;  // 节点力向量 (瞬时存储)
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> solution_rf; // 约束力向量

		// --- 系统相关数据 ---
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> inv_mass_matrix; // 节点质量矩阵逆 (对角矩阵)
		Material_elastic mMatElastic;                    // 材料属性 (目前只有一个)

		// --- FEM 计算辅助数据 ---
		Gauss_integral_1d gauss_intg_1d_1;              // 高斯积分值 (1D, 1次)
		FEValues<8, 3> feValuesHex;                     // FE 值类 (所有单元, 节点自由度)

		// --- 时间相关数据 ---
		Types::Real time;                   // 前一步模型时间
		Types::Real totalTime;              // 模型时间
		Types::Real dt_i;                   // 前一步时间步长 (dt_n+1/2)
		Types::Real dt_i_1;                 // 上一时间步长 (dt_n-1/2)
		Types::Real factor_timeStep;        // 时间步长比例系数

		// --- 模型参数数据 ---
		Types::Real c_dr;                   // 沙漏阻尼系数
		Types::Real c_cr;                   // 沙漏粘性系数
		Types::Real Cvisl;                  // 体积粘性系数
		Types::Real Cvisq;                  // 体积粘性稳定系数

		// <<< 新增：体积变化率计算阈值 >>>
		const Types::Real volRateThreshold = static_cast<Types::Real>(1e-10);

		// --- 重力、边界条件初始化数据 ---
		std::tuple<bool, std::string, Types::Real, Types::Real, Types::Real, Types::Real> gravity; // 重力
		std::vector<std::tuple<std::string, std::string, std::string, Types::Real>> dsload;        // 分布外力 (瞬时存储)
		std::vector<std::pair<std::vector<std::size_t>, std::array<std::size_t, 3>>> boundary_spc_node; // 位移约束节点集合
		std::vector<std::pair<std::vector<std::size_t>, std::pair<std::array<std::size_t, 2>, Types::Real>>> boundary_vel_node; // 速度约束节点集合
		std::vector<std::pair<std::vector<std::size_t>, std::pair<std::size_t, Types::Real>>> ini_vel_generation; // 初始速度生成节点集合
		std::map<int, int> mPart_PID_MID;              // Part ID 到 Material ID 映射 (瞬时存储)
		std::map<std::string, std::vector<std::size_t>> map_set_node_list; // 节点集合到节点数组映射
		std::map<unsigned int, unsigned int> map_set_part_ID; // Part Set ID 映射? (瞬时存储)
		std::map<std::string, std::vector<std::array<Types::Real, 2>>> map_amp_list; // 幅值曲线生成参数数组到数据数组映射
		std::map<std::string, std::vector<std::size_t>> map_set_ele_list; // 单元集合到单元数组映射
		std::map<std::string, std::pair<std::string, std::string>> map_set_surface_list; // 表面集合 (单元, 表面编号) 映射? (瞬时存储)

		// --- 输出相关数据 ---
		Types::Real time_interval;          // 时间间隔
		Types::Real time_output;            // 上一时间输出时间
		bool time_interval_set;             // 标记time_interval是否已设置

		// --- 计算辅助数据 (单元) ---
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> von_mises; // Von Mises 应力
		Eigen::Matrix<Types::Real, Eigen::Dynamic, 1> dtStable;   // 每单元粘性/声速稳定时间步

		// --- FSI 相关数据 ---
		// (存储 FSI 边界信息，用于构建结构体元素映射)
		std::vector<std::pair<int, std::vector<int>>> fsiBoundaryElementFaces; // FSI 边界单元面 -> 边界`面数组
		std::map<int, std::vector<std::array<Types::Real, 3>>> map_eleIndex_virParticle_unit_data; // 单元 -> 虚拟粒子单位坐标数组
		std::map<int, std::vector<int>> map_eleIndex_virParticles_index; // 单元 -> 虚拟粒子全索引数组

		// --- 新增变量 ---
		std::size_t currentStepIndex; // 新增变量：当前时间步骤索引
		std::vector<StepData> steps;  // 新增变量：所有步骤的数据

	};

} // namespace EnSC