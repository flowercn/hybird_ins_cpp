cmake_minimum_required(VERSION 3.10)

# 查找 Eigen3 (必须)
find_package(Eigen3 REQUIRED)

# 定义源文件
set(PSINS_SRC
    glv.cpp
    earth.cpp
    ins_math.cpp
    ins_state.cpp
    kf_state.cpp
    sins_engine.h
)

# 生成静态库 psins_core
add_library(psins_core STATIC ${PSINS_SRC})

# 配置头文件搜索路径
# PUBLIC 属性确保链接了这个库的目标也能自动找到这些头文件
target_include_directories(psins_core PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    ${EIGEN3_INCLUDE_DIR}
)

# 链接 Eigen3
target_link_libraries(psins_core PUBLIC Eigen3::Eigen)