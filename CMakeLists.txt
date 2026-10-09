cmake_minimum_required(VERSION 3.15...3.27)

# Prefer a compiler that links against the system's own glibc over conda's own. Activating a conda
# environment that has gxx_linux-64/gcc_linux-64 installed auto-sets CXX/CC to conda's own compiler
# (e.g. .../envs/<env>/bin/x86_64-conda-linux-gnu-c++) -- it supports C++20, but links against
# conda's own bundled sysroot package, which can target a newer glibc than the machine running the
# code actually has (this has bitten us: the extension built fine, then failed at import with
# "GLIBC_2.32 not found" on a cluster login node still on glibc 2.28). RHEL's gcc-toolset packages
# give a newer GCC while still linking against the real system glibc, so prefer the newest one
# available whenever CXX/CC is unset or points inside the active conda environment -- i.e. whenever
# it looks like conda's own auto-set compiler rather than something the caller deliberately chose
# (an explicit -DCMAKE_CXX_COMPILER=..., or a CXX pointing outside the conda environment, such as a
# manually module-loaded toolset, is always respected as-is). On a machine without
# /opt/rh/gcc-toolset-* (e.g. not RHEL), this is a no-op and CMake's normal discovery proceeds.
if(DEFINED ENV{CONDA_PREFIX} AND NOT DEFINED CMAKE_CXX_COMPILER)
  set(_cxx_is_conda_or_unset TRUE)
  if(DEFINED ENV{CXX} AND NOT "$ENV{CXX}" MATCHES "^$ENV{CONDA_PREFIX}/")
    set(_cxx_is_conda_or_unset FALSE)
  endif()
  if(_cxx_is_conda_or_unset)
    file(GLOB _gcc_toolsets "/opt/rh/gcc-toolset-*")
    if(_gcc_toolsets)
      list(SORT _gcc_toolsets COMPARE NATURAL ORDER DESCENDING)
      list(GET _gcc_toolsets 0 _latest_toolset)
      if(EXISTS "${_latest_toolset}/root/usr/bin/g++" AND EXISTS "${_latest_toolset}/root/usr/bin/gcc")
        message(STATUS "Building inside a conda environment: using ${_latest_toolset} instead of "
                        "conda's own compiler, so the extension links against the system's glibc "
                        "rather than conda's bundled (newer) one.")
        set(CMAKE_C_COMPILER "${_latest_toolset}/root/usr/bin/gcc")
        set(CMAKE_CXX_COMPILER "${_latest_toolset}/root/usr/bin/g++")
      endif()
    endif()
    unset(_gcc_toolsets)
    unset(_latest_toolset)
  endif()
  unset(_cxx_is_conda_or_unset)
endif()

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g -O0 -ggdb3)
  set(CMAKE_STRIP_COMMAND "")
  set(CMAKE_STRIP_FLAG_DEBUG "")
  set(CMAKE_STRIP_FLAG_RELEASE "")
  add_compile_definitions(DEBUG_CHECKS)
endif()

set(CMAKE_CXX_STANDARD 20)

pybind11_add_module(_atmosphere_bgs MODULE src/main.cpp src/laguerre_diagram.cpp src/discretized_line_segment.cpp src/rasterizer.cpp)
target_include_directories(_atmosphere_bgs PUBLIC extern/eigen extern/pybind11/include extern/fmt/include extern/cpp-TimSort/include)

install(TARGETS _atmosphere_bgs LIBRARY DESTINATION .)
