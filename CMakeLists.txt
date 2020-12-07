cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(ttial)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)
enable_testing()

include(FetchContent)

## GCC
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code -Wstrict-aliasing -O2)
## Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wunreachable-code -Wstrict-aliasing -O2)
endif()

# download Eigen
FetchContent_Declare(
    eigen
    URL "https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz"
    URL_HASH MD5=9e30f67e8531477de4117506fe44669b)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
endif()

include_directories(SYSTEM ${eigen_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

add_executable(ttial ttial.cpp)
