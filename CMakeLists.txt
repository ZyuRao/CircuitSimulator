cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(CircuitSimulator LANGUAGES CXX)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 设置 Eigen 库的路径（替换为你的安装路径）
set(EIGEN3_INCLUDE_DIRS "C:/mingw64/include/eigen3")

# 将 include 文件夹加入头文件搜索路径
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIRS}    
)

# 添加源代码文件
file(GLOB SRC_FILES
    src/*.cpp
)
# 添加可执行文件
add_executable(${PROJECT_NAME} ${SRC_FILES})

