#pragma once
#include <string>
#include <vector>
#include <map>
#include <array>
#include "types.h"

namespace EnSC {

// 单个Part结构
struct InpPart {
    std::string name; // Part名称
    std::vector<Types::Point<3>> nodes; // 节点坐标
    std::vector<std::vector<Types::Vertex_index>> elements; // 单元（节点索引）
    std::map<std::string, std::vector<std::size_t>> node_sets; // 节点集
    std::map<std::string, std::vector<std::size_t>> element_sets; // 单元集
    std::string material_name; // 材料名
    // 可扩展：表面、section等
};

// Instance结构（装配中的一个部件实例）
struct InpInstance {
    std::string name;      // 实例名
    std::string part_name; // 引用的Part名
    // 新增：坐标变换相关
    std::array<Types::Real, 3> translation; // 平移向量
    // 新增：节点和单元索引偏移
    std::size_t node_start_index;    // 全局节点起始索引
    std::size_t element_start_index; // 全局单元起始索引
};

// Assembly结构
struct InpAssembly {
    std::string name; // 装配名
    std::vector<InpInstance> instances; // 实例列表
    // 修改：节点集和单元集的存储结构
    // 从 map<string, vector<size_t>> 改为 map<string, map<string, vector<size_t>>>
    // 第一个string是集合名，第二个string是实例名，vector是该实例中的索引列表
    std::map<std::string, std::map<std::string, std::vector<std::size_t>>> node_sets; // 节点集（装配级）
    std::map<std::string, std::map<std::string, std::vector<std::size_t>>> element_sets; // 单元集（装配级）
    // 可扩展：表面、section等
};

// 材料结构
struct InpMaterial {
    std::string name;
    double density;
    double E;
    double v;
    // 可扩展：更多材料参数
};

// 顶层InpData结构
struct InpData {
    std::vector<InpPart> parts;
    std::vector<InpAssembly> assemblies;
    std::vector<InpMaterial> materials;
    // 可扩展：step、output等
};

} // namespace EnSC 