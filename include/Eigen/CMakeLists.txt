cmake_minimum_required(VERSION 3.10)
project(hybrid_ins_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 REQUIRED)
find_package(GTSAM REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${GTSAM_INCLUDE_DIR})
include_directories(include)

# 纯惯导与融合引擎 (去掉了已删除的 ins_align.cpp)
add_executable(fgo_engine src/fgo_main.cpp src/ins_core.cpp)
target_link_libraries(fgo_engine gtsam)

# 生成 CAI 数据工具
add_executable(gen_cai src/gen_cai_data.cpp src/ins_core.cpp)
