#include <gtest/gtest.h>
#include "../include/exDyna3D.h"

namespace EnSC {
namespace testing {

// 创建一个继承自exDyna3D的测试类，以访问protected成员
class TestableExDyna3D : public exDyna3D {
public:
    // 公开apply_fsiSph_nodeForce方法供测试使用
    using exDyna3D::apply_fsiSph_nodeForce;
    
    // 公开需要访问的protected成员
    using exDyna3D::vertices;
    using exDyna3D::hexahedron_elements;
    using exDyna3D::solution_a;
    using exDyna3D::map_eleIndex_virParticle_unit_data;
    using exDyna3D::map_eleIndex_virParticles_index;
    using exDyna3D::fsi_share_data;
};

// 创建一个测试fixture用于测试apply_fsiSph_nodeForce函数
class ApplyFsiSphNodeForceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 设置测试用例的初始状态
        setupTestData();
    }

    void TearDown() override {
        // 测试结束后的清理工作
    }

    void setupTestData() {
        // 创建一个简化的测试模型
        // 1. 设置一个单元的网格
        exdyna.hexahedron_elements.resize(1);  // 一个六面体单元
        exdyna.vertices.resize(8);  // 8个节点

        // 设置节点坐标（简单的单位立方体）
        exdyna.vertices[0] = {0.0, 0.0, 0.0};
        exdyna.vertices[1] = {1.0, 0.0, 0.0};
        exdyna.vertices[2] = {1.0, 1.0, 0.0};
        exdyna.vertices[3] = {0.0, 1.0, 0.0};
        exdyna.vertices[4] = {0.0, 0.0, 1.0};
        exdyna.vertices[5] = {1.0, 0.0, 1.0};
        exdyna.vertices[6] = {1.0, 1.0, 1.0};
        exdyna.vertices[7] = {0.0, 1.0, 1.0};

        // 设置单元顶点索引
        std::vector<Types::Vertex_index> vertexIndices = {0, 1, 2, 3, 4, 5, 6, 7};
        exdyna.hexahedron_elements[0].set_verticesIndex(vertexIndices);

        // 2. 创建虚拟粒子数据
        // 添加一个虚拟粒子的单位坐标（在参考单元内）
        std::array<Types::Real, 3> unitCoord = {0.0, 0.0, 0.0};  // 例如在单元中心
        exdyna.map_eleIndex_virParticle_unit_data[0].push_back(unitCoord);

        // 添加虚拟粒子索引映射
        exdyna.map_eleIndex_virParticles_index[0].push_back(0);

        // 3. 设置虚拟粒子力数据
        exdyna.fsi_share_data.FSI_virtualParticles_nodeForce.resize(1);
        exdyna.fsi_share_data.FSI_virtualParticles_nodeForce[0] = {10.0, 20.0, 30.0};  // x, y, z方向力

        // 4. 初始化解向量
        exdyna.solution_a = Eigen::Matrix<Types::Real, Eigen::Dynamic, 1>::Zero(24);  // 8个节点，每个节点3个自由度
    }

    // 创建一个测试用的TestableExDyna3D实例
    TestableExDyna3D exdyna;
};

// 测试apply_fsiSph_nodeForce函数的主要功能
TEST_F(ApplyFsiSphNodeForceTest, AppliesVirtualParticleForces) {
    // 调用被测试的函数
    exdyna.apply_fsiSph_nodeForce();

    // 验证力已正确应用到节点
    // 检查所有节点的加速度向量是否包含了虚拟粒子的力贡献
    for (int i = 0; i < 8; i++) {
        int dof0 = 3 * i;
        // 由于形函数的值在(0,0,0)对所有节点的贡献通常是均匀的，所以期望所有节点受到均匀的力
        // 具体值取决于六面体单元的形函数实现
        EXPECT_GT(exdyna.solution_a[dof0], 0.0) << "X方向力应该被应用到节点 " << i;
        EXPECT_GT(exdyna.solution_a[dof0 + 1], 0.0) << "Y方向力应该被应用到节点 " << i;
        EXPECT_GT(exdyna.solution_a[dof0 + 2], 0.0) << "Z方向力应该被应用到节点 " << i;
    }
}

// 测试边界条件
TEST_F(ApplyFsiSphNodeForceTest, HandlesEmptyVirtualParticles) {
    // 清空虚拟粒子数据
    exdyna.map_eleIndex_virParticle_unit_data.clear();
    exdyna.map_eleIndex_virParticles_index.clear();
    
    // 将解向量置零
    exdyna.solution_a.setZero();
    
    // 调用函数
    exdyna.apply_fsiSph_nodeForce();
    
    // 验证没有力被应用
    for (int i = 0; i < exdyna.solution_a.size(); i++) {
        EXPECT_DOUBLE_EQ(exdyna.solution_a[i], 0.0) << "没有虚拟粒子时不应该应用力";
    }
}

// 测试中心粒子受到单位力沿着x方向的情况
TEST_F(ApplyFsiSphNodeForceTest, AppliesUnitForceInXDirection) {
    // 重置加速度向量
    exdyna.solution_a.setZero();
    
    // 修改虚拟粒子力为沿x方向的单位力
    exdyna.fsi_share_data.FSI_virtualParticles_nodeForce[0] = {1.0, 0.0, 0.0};  // 只在x方向有单位力
    
    // 调用被测试的函数
    exdyna.apply_fsiSph_nodeForce();
    
    // 打印所有节点的x方向力，用于调试
    std::cout << "节点的x方向力分布：" << std::endl;
    for (int i = 0; i < 8; i++) {
        int dof0 = 3 * i;
        std::cout << "节点 " << i << ": " << exdyna.solution_a[dof0] << std::endl;
    }
    
    // 验证力已正确应用到节点
    for (int i = 0; i < 8; i++) {
        int dof0 = 3 * i;
        // 检查x方向应该有力
        EXPECT_GT(exdyna.solution_a[dof0], 0.0) << "X方向力应该被应用到节点 " << i;
        // y和z方向应该没有力
        EXPECT_DOUBLE_EQ(exdyna.solution_a[dof0 + 1], 0.0) << "Y方向不应该有力在节点 " << i;
        EXPECT_DOUBLE_EQ(exdyna.solution_a[dof0 + 2], 0.0) << "Z方向不应该有力在节点 " << i;
    }
    
    // 验证力的总和等于输入的力（守恒）
    Types::Real totalForceX = 0.0;
    for (int i = 0; i < 8; i++) {
        totalForceX += exdyna.solution_a[3 * i];
    }
    EXPECT_NEAR(totalForceX, 1.0, 1e-10) << "X方向的总力应该等于输入的单位力";
}

}  // namespace testing
}  // namespace EnSC 