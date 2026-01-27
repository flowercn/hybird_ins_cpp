cmake_minimum_required(VERSION 3.10)

# 查找 Eigen3
find_package(Eigen3 REQUIRED)

# 定义源文件
# [修正]：将 sins_engine.h 替换为 sins_engine.cpp
set(PSINS_SRC
    earth.cpp
    ins_math.cpp
    ins_state.cpp
    kf_state.cpp
    sins_engine.cpp 
)

# 生成静态库 psins_core
add_library(psins_core STATIC ${PSINS_SRC})

# 配置头文件搜索路径
target_include_directories(psins_core PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    ${EIGEN3_INCLUDE_DIR}
)

# 链接 Eigen3
target_link_libraries(psins_core PUBLIC Eigen3::Eigen)