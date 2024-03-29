cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg --coverage")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(startle VERSION 1.0)

add_subdirectory(src/)

target_include_directories(startle PUBLIC "${startle_SOURCE_DIR}/include")
target_include_directories(startle PUBLIC "${startle_SOURCE_DIR}/third-party/csv-parser/include")

add_subdirectory(third-party/spdlog)
add_subdirectory(third-party/pprint)
add_subdirectory(third-party/csv-parser)
add_subdirectory(third-party/argparse)
add_subdirectory(third-party/json)

enable_testing()
