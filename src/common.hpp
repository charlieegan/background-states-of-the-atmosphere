#ifndef COMMON_HPP
#define COMMON_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
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
#include <queue>
#include <deque>

#include <exception>
#include <functional>
#include <optional>

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

[[maybe_unused]] static std::partial_ordering operator*(const std::partial_ordering &o0, const std::partial_ordering &o1) {
  if (o0 == 0 || o1 == 0)
    return 0 <=> 0;
  if ((o0 < 0 && o1 < 0) || (o0 > 0 && o1 > 0))
    return 1 <=> 0;
  return 0 <=> 1;
}

[[maybe_unused]] static std::string to_string(const std::partial_ordering &o) {
  if (o == std::partial_ordering::less)
    return "less";
  if (o == std::partial_ordering::equivalent)
    return "equivalent";
  if (o == std::partial_ordering::greater)
    return "greater";
  return "unordered";
}

template <typename T>
struct type_name {
  static std::string value() {
    if (std::is_floating_point_v<T>)
      return FORMAT("float{}", 8 * sizeof(T));
    else if (std::is_integral_v<T> && std::is_signed_v<T>)
      return FORMAT("int{}", 8 * sizeof(T));
    else if (std::is_integral_v<T> && std::is_unsigned_v<T>)
      return FORMAT("uint{}", 8 * sizeof(T));
    return typeid(T).name();
  }
};

template <typename T>
[[maybe_unused]] bool significantly_less(const T &lhs, const T &rhs,
                                        const T &rel = 100 * std::numeric_limits<T>::epsilon(),
                                        const T &abs = 100 * std::numeric_limits<T>::epsilon()) {
  T mean = 0.5 * (lhs + rhs);
  T eps = std::max(abs, std::abs(mean) * rel);
  return lhs < rhs - eps;
}

#endif
