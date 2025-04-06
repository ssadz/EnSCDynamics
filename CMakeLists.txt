# CMake 最低版本要求
cmake_minimum_required(VERSION 3.12)

# 这个项目是一个基于EnSC框架的显式动力学求解器，使用Eigen进行数值计算，VTK进行可视化，通过OpenMP实现并行加速
# 项目名称和语言
project(EnSCDynamics LANGUAGES CXX)

# --- 编译器设置 (优先在 Windows 上使用 Clang) ---
if(WIN32)
    message(STATUS "检测到 Windows 平台，尝试查找并使用 Clang...")
    find_program(CLANG_C_EXECUTABLE NAMES clang clang.exe DOC "Clang C Compiler")
    find_program(CLANG_CXX_EXECUTABLE NAMES clang++ clang++.exe DOC "Clang C++ Compiler")

    if(CLANG_C_EXECUTABLE AND CLANG_CXX_EXECUTABLE)
        # 设置 Clang 为编译器
        set(CMAKE_C_COMPILER "${CLANG_C_EXECUTABLE}")
        set(CMAKE_CXX_COMPILER "${CLANG_CXX_EXECUTABLE}")
        message(STATUS "找到并设置 Clang 作为编译器:")
        message(STATUS "  C Compiler: ${CMAKE_C_COMPILER}")
        message(STATUS "  CXX Compiler: ${CMAKE_CXX_COMPILER}")

        # --- 添加 Clang 路径到 CMAKE_PREFIX_PATH ---
        # 这对于查找与 Clang 安装相关的库（如其自带的 OpenMP 运行时）可能很重要
        get_filename_component(CLANG_BIN_DIR "${CLANG_C_EXECUTABLE}" DIRECTORY)
        get_filename_component(CLANG_ROOT_DIR "${CLANG_BIN_DIR}" DIRECTORY) # 获取 bin 目录的上级目录
        # 对于标准的 LLVM 安装，这通常是正确的根目录
        # 如果你的 Clang 安装结构不同，可能需要调整 get_filename_component 的层级

        if(IS_DIRECTORY "${CLANG_ROOT_DIR}")
            list(APPEND CMAKE_PREFIX_PATH "${CLANG_ROOT_DIR}")
            message(STATUS "将 Clang 根目录 ${CLANG_ROOT_DIR} 添加到 CMAKE_PREFIX_PATH 以帮助查找相关库 (如 OpenMP)。")
            message(STATUS "更新后的 CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
        else()
             message(WARNING "未能确定有效的 Clang 根目录 (${CLANG_ROOT_DIR})，可能影响查找 Clang 相关库。")
        endif()
        # --- Clang 路径添加结束 ---

    else()
        message(WARNING "未完整找到 Clang (clang 和 clang++)。将使用 CMake 自动检测到的其他编译器 (如 MSVC)。")
    endif()
else()
    message(STATUS "非 Windows 平台，使用系统默认 C/C++ 编译器。")
endif()

# --- 项目通用设置 ---
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "构建类型 (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
    message(STATUS "未指定构建类型，已自动设置为: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "当前构建类型: ${CMAKE_BUILD_TYPE}")
endif()

# --- Release模式下的编译优化选项 ---
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # GCC或Clang编译器
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
    message(STATUS "为Release模式设置GCC/Clang优化选项: -O3 -DNDEBUG")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # MSVC编译器
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG /arch:AVX2 /GL")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    message(STATUS "为Release模式设置MSVC优化选项: /O2 /DNDEBUG /arch:AVX2 /GL /LTCG")
endif()

# --- 添加本地依赖库的搜索路径 (继续追加) ---
# VTK 路径
if(WIN32)
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/VTK_win64")
    set(LOCAL_VTK_BIN_DIR "${CMAKE_SOURCE_DIR}/VTK_win64/bin")
    message(STATUS "已添加 Windows VTK 路径到搜索路径: ${CMAKE_SOURCE_DIR}/VTK_win64")
else()
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/VTK")
    set(LOCAL_VTK_BIN_DIR "")
    message(STATUS "已添加 Unix-like VTK 路径到搜索路径: ${CMAKE_SOURCE_DIR}/VTK")
endif()

# Eigen 路径
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/eigen")
message(STATUS "已添加 Eigen 路径到搜索路径: ${CMAKE_SOURCE_DIR}/eigen")

message(STATUS "最终的库搜索路径 (CMAKE_PREFIX_PATH): ${CMAKE_PREFIX_PATH}")

# 添加 Explicit3D 子目录
add_subdirectory(Explicit3D)

# --- 构建可执行文件并链接 Explicit3D 对象库 ---
add_executable(${PROJECT_NAME} main.cpp)

# 为主程序添加包含路径
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    Explicit3D
)

# --- 构建后操作 (仅 Windows: 复制 VTK DLLs) ---
if(WIN32 AND LOCAL_VTK_BIN_DIR)
    if(EXISTS "${LOCAL_VTK_BIN_DIR}")
        message(STATUS "配置构建后命令：将 VTK DLLs 从 ${LOCAL_VTK_BIN_DIR} 复制到输出目录...")
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
                    "${LOCAL_VTK_BIN_DIR}/"
                    "$<TARGET_FILE_DIR:${PROJECT_NAME}>/"
            COMMENT "正在将 VTK DLLs 从 ${LOCAL_VTK_BIN_DIR} 复制到 $<TARGET_FILE_DIR:${PROJECT_NAME}>"
            VERBATIM
        )
    else()
        message(WARNING "指定的 VTK Bin 目录 (${LOCAL_VTK_BIN_DIR}) 不存在，无法配置 DLL 复制。请检查路径是否正确。")
    endif()
endif()

# --- 为所有平台复制项目配置文件和数据目录 ---
# 复制 project.txt 到输出目录
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/project.txt"
            "$<TARGET_FILE_DIR:${PROJECT_NAME}>/project.txt"
    COMMENT "正在复制 project.txt 到 $<TARGET_FILE_DIR:${PROJECT_NAME}>"
    VERBATIM
)

# 复制 inp 目录到输出目录
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_FILE_DIR:${PROJECT_NAME}>/inp"
    COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
            "${CMAKE_SOURCE_DIR}/inp"
            "$<TARGET_FILE_DIR:${PROJECT_NAME}>/inp"
    COMMENT "正在复制 inp 目录到 $<TARGET_FILE_DIR:${PROJECT_NAME}>"
    VERBATIM
)

# --- Release构建说明 ---
# 要进行Release构建，请使用以下命令：
#
# Windows上（命令行）:
# cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
# cmake --build build --config Release
#
# Windows上（Visual Studio）:
# 在生成解决方案后，将构建配置从"Debug"切换为"Release"
#
# Linux/Mac上:
# cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
# cmake --build build
#
# 这将启用所有优化选项，包括-O3、-DNDEBUG和针对本机的CPU指令集优化