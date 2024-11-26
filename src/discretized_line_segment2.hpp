#ifndef DISCRETIZED_LINE_SEGMENT
#define DISCRETIZED_LINE_SEGMENT

#include "common.hpp"
#include "physical_parameters.hpp"
#include "simulation_parameters.hpp"

template <typename T>
class discretized_line_segment
{
public:
  typedef Eigen::Matrix<T, Eigen::Dynamic, 2, Eigen::RowMajor> Matrix2X;
  typedef Eigen::Vector2<T> Vector2;
  typedef Eigen::VectorX<T> VectorX;

  const physical_parameters *phys;
  int max_resolution;
  Vector2 aspect;

  // start, end and direction in linear coords
  Vector2 start, end, direction;
  
  // the coordinate lists are all on the same discretizaion
  // s is one shorter than the rest, since there is no intersect after the last point
  VectorX lams; // positions in [0,1] at which discretization points are (linearly spaced in linear domain)
  Matrix2X x; // coordinates in s-p-coords
  Matrix2X t; // tangent vectors in s-p-coords
  Matrix2X s; // intersections between current tangent with next tangent in s-p-coords

  T errb; // area error bound
  T area; // (approx) area below (discretized) line segment
  
  discretized_line_segment(const Eigen::Ref<const Vector2> &zeta_s, // start point in linear coords
                           const Eigen::Ref<const Vector2> &zeta_e, // end point in linear coords
                           const physical_parameters &phys,
                           const simulation_parameters &sim);

  Vector2 get_z(const T &lam) const; // get linear coord at position lam in [0,1] along line
  Vector2 get_x(const T &lam) const; // get s-p-coords at position lam in [0,1] along line
  Vector2 get_t(const T &lam) const; // get s-p-tangent at position lam in [0,1] along line

  void refine(int newlen = -1); // refine resolution to newlen (if newlen < 2, roughly double it)

  int size() const; // number of points used for discretization
  
  std::string repr() const;

  static void bind(py::module_ &m);
};

#endif
