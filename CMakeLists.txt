cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(CircuitSimulator)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)

# 设置 Eigen 库的路径（替换为你的安装路径）
set(Eigen_3DIR "D:/libraries/eigen-3.4.0/cmake")

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
# 将 include 文件夹加入头文件搜索路径
include_directories(${PROJECT_SOURCE_DIR}/include)

# 添加源代码文件
set(SOURCE_FILES
    src/circuit.cpp
    src/dcanalysis.cpp
    src/main.cpp
    src/mna.cpp
    src/solver.cpp
    src/sstateanalysis.cpp
    src/tanalisis.cpp
    src/utils.cpp
)

# 添加可执行文件
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
