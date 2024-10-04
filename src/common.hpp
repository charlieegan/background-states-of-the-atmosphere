#ifndef COMMON_HPP
#define COMMON_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
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
#include <set>
#include <vector>

#include <exception>
#include <functional>

#include <random>

#include <concepts>

#include "timer.hpp"

#if __has_include(<format>)
#include <format>
#define FORMAT std::format
#else
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#define FORMAT fmt::format
#endif

static inline auto sqr(const auto &x) { return x * x; }

static std::partial_ordering operator*(const std::partial_ordering &o0, const std::partial_ordering &o1) {
  if (o0 == 0 || o1 == 0)
    return 0 <=> 0;
  if ((o0 < 0 && o1 < 0) || (o0 > 0 && o1 > 0))
    return 1 <=> 0;
  return 0 <=> 1;
}

static std::string to_string(const std::partial_ordering &o) {
  if (o == std::partial_ordering::less)
    return "less";
  if (o == std::partial_ordering::equivalent)
    return "equivalent";
  if (o == std::partial_ordering::greater)
    return "greater";
  return "unordered";
}

#endif
