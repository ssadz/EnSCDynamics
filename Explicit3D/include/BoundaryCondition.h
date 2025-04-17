#pragma once

#include <vector>
#include <array>
#include <utility> // For std::pair
#include <cstddef> // For std::size_t
#include "types.h" // For Types::Real

namespace EnSC {

// --- 边界条件结构体 ---
struct BoundaryCondition {
    // 位移约束节点集合
    std::vector<std::pair<std::vector<std::size_t>, std::array<std::size_t, 3>>> spc_nodes;
    
    // 速度约束节点集合
    std::vector<std::pair<std::vector<std::size_t>, std::pair<std::array<std::size_t, 2>, Types::Real>>> vel_nodes;
};

} // namespace EnSC 