cmake_minimum_required(VERSION 3.12)
project(ExampleProject)

# Find the Python3 package
find_package(Python3 REQUIRED)
find_package (Python3 COMPONENTS Interpreter Development)

include_directories(include)

# Add your C++ executable
add_executable(example example.cpp)

# Include Python3 headers
target_include_directories(example PRIVATE ${Python3_INCLUDE_DIRS})

# Link against Python3 library
target_link_libraries(example PRIVATE Python3::Python)