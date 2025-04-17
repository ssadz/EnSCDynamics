#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <utility> // for std::get

#include "types.h"
#include "BoundaryCondition.h"

namespace EnSC {

// --- 步骤数据结构 ---
struct StepData {
    std::string name;                 // 步骤名称
    Types::Real timePeriod;           // 时间周期
    BoundaryCondition boundary;       // 边界条件
    bool resetSpcBoundary;            // 是否重置了位移边界条件 (Boundary, op=NEW)
    bool resetVelBoundary;            // 是否重置了速度边界条件 (Boundary, op=NEW, type=VELOCITY)
    
    // 载荷相关
    std::tuple<bool, std::string, Types::Real, Types::Real, Types::Real, Types::Real> gravity; // 重力
    std::vector<std::tuple<std::string, std::string, std::string, Types::Real>> dsload;        // 分布外力
    bool resetDload;                  // 是否重置了分布载荷 (DLOAD, op=NEW)
    bool resetDsload;                 // 是否重置了分布面载荷 (DSLOAD, op=NEW)
    
    // 构造函数，初始化默认值
    StepData() : timePeriod(0.0), resetSpcBoundary(false), resetVelBoundary(false), 
                 resetDload(false), resetDsload(false) {
        // 初始化gravity元组的第一个元素（使能标志）为false
        std::get<0>(gravity) = false;
    }
    // 其他特定于步骤的数据可以在这里添加
};

} // namespace EnSC 