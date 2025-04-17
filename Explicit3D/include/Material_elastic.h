#pragma once

#include "types.h"

namespace EnSC {

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

} // namespace EnSC 