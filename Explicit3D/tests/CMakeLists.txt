# 测试目录的 CMake 配置

# 收集所有测试源文件
file(GLOB TEST_SOURCES "*.cpp")

# 创建测试可执行文件
add_executable(Explicit3DTests ${TEST_SOURCES})

# 设置测试可执行文件的输出目录，与主项目相同
set_target_properties(Explicit3DTests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}"
)

# 链接到 Explicit3D 库和 GoogleTest
target_link_libraries(Explicit3DTests PRIVATE
    Explicit3D
    gtest
    gtest_main
    VTK::CommonCore
    VTK::CommonDataModel
    VTK::IOXML
    Eigen3::Eigen
)

# 为测试可执行文件添加包含目录
target_include_directories(Explicit3DTests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# 添加测试到 CTest
add_test(NAME Explicit3DTests COMMAND Explicit3DTests) 