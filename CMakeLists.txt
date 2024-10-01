cmake_minimum_required(VERSION 3.15...3.27)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

set (CMAKE_CXX_STANDARD 20)

include_directories(extern/eigen extern/pybind11/include extern/fmt/include extern/cpp-TimSort/include)

pybind11_add_module(dvector SHARED src/dvector.cpp)

pybind11_add_module(_atmosphere_bgs src/main.cpp)
target_link_libraries(_atmosphere_bgs PUBLIC dvector)

install(TARGETS _atmosphere_bgs LIBRARY DESTINATION .)
