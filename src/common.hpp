#ifndef COMMON_HPP
#define COMMON_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <Eigen/Dense>

#include <cmath>
#include <algorithm>
#include <numeric>

#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

#include <map>
#include <list>
#include <vector>

#include <exception>
#include <functional>

#include <random>

#include "timer.hpp"

#if __has_include(<format>)
#include <format>
#define FORMAT std::format
#else
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#define FORMAT fmt::format
#endif



#endif
