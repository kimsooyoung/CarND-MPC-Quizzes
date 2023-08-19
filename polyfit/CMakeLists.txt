project(q1)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

set(sources src/main.cpp)

find_package(Eigen3 REQUIRED)
find_package(Python3 REQUIRED)
find_package (Python3 COMPONENTS Interpreter Development)

# include_directories(/usr/local/include)
# include_directories(src/Eigen-3.3)
include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
)

link_directories(/usr/local/lib)

add_executable(third_poly src/third_poly.cpp)
add_executable(polyfit src/main.cpp)

target_link_libraries(third_poly Python3::Python)
target_link_libraries(polyfit Python3::Python)
