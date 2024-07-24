#ifndef TRANSFORM2D_HPP
#define TRANSFORM2D_HPP

#include <Eigen/Dense>

class transform2d {

  // differential transform: transform point x and tangent vector v to point and tangent vector in output coordinates
  virtual std::pair<Eigen::Vector2d, Eigen::Vector2d>
  dtf(const Eigen::Ref<const Eigen::Vector2d> &x, const Eigen::Ref<const Eigen::Vector2d> &v);
}

#endif
