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
#include <forward_list>
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

#if __has_include(<stacktrace>) and __cplusplus >= 202302L
#include <stacktrace>
std::string get_trace() {
  return std::to_string(std::stacktrace::current());
}
#else
std::string get_trace() {
  return "use C++23 for stacktrace support";
}
#endif

inline auto sqr(const auto &x) { return x * x; }

std::partial_ordering operator*(const std::partial_ordering &o0, const std::partial_ordering &o1) {
  if (o0 == 0 || o1 == 0)
    return 0 <=> 0;
  if ((o0 < 0 && o1 < 0) || (o0 > 0 && o1 > 0))
    return 1 <=> 0;
  return 0 <=> 1;
}

std::string to_string(const std::partial_ordering &o) {
  if (o == std::partial_ordering::less)
    return "less";
  if (o == std::partial_ordering::equivalent)
    return "equivalent";
  if (o == std::partial_ordering::greater)
    return "greater";
  return "unordered";
}

#endif
