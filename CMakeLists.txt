cmake_minimum_required(VERSION 3.15...3.27)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)

pybind11_add_module(_atmosphere_bgs MODULE src/main.cpp src/laguerre_diagram.cpp src/discretized_line_segment2.cpp src/rasterizer.cpp)
target_include_directories(_atmosphere_bgs PUBLIC extern/eigen extern/pybind11/include extern/fmt/include extern/cpp-TimSort/include)

install(TARGETS _atmosphere_bgs LIBRARY DESTINATION .)
