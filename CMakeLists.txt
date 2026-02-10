cmake_minimum_required(VERSION 3.10)
project(Sandym)

set(CMAKE_C_STANDARD 99)


# 添加可执行文件，列出所有源文件
add_executable(Sandparticlesimulator main.c SandSim.c)
target_link_libraries(Sandparticlesimulator m)
