cmake_minimum_required(VERSION 3.12)
project(WASCO CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization and warning flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -O2)

# Common source files
set(COMMON_SOURCES
    graph.cpp
    functions.cpp
    exp_func.cpp
    experiment.cpp
    exact.cpp
    compare.cpp
)

# Main wasco binary
add_executable(wasco
    main.cpp
    ${COMMON_SOURCES}
)

# Include current directory for headers
target_include_directories(wasco PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
