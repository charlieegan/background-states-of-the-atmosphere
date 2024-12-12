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

  // start, end and direction in linear coords
  Vector2 start; //!< start of segment in linear coords
  Vector2 end; //!< end of segment in linear coords
  Vector2 direction; //!< normalized direction of segment in linear coords

  std::shared_ptr<physical_parameters> phys; //!< physical parameters defining the coordante transform
  int max_resolution; //!< maximal resolution (number of points) for the approximation
  Vector2 aspect; //!< coordinate stretch to apply before norm for even spacing (to correct for aspect ratio far off 1)
  
  // the coordinate lists are all on the same discretizaion
  // s is one shorter than the rest, since there is no intersect after the last point
  VectorX lams; //!< positions in [0,1] at which discretization points are (linearly spaced in linear domain)
  Matrix2X x; //!< coordinates in s-p-coords
  Matrix2X t; //!< tangent vectors in s-p-coords
  Matrix2X s; //!< intersections between current tangent with next tangent in s-p-coords

  T errb; //!< area error bound
  T area; //!< (approx) area below (discretized) line segment in s-p-coords

  /*! construct discretized_line_segment between two points given in linear coords
   *
   * \param zeta_s start point in linear coords
   * \param zeta_e end point in linear coords
   * \param phys physical parameters defining transform
   * \param sim simulation parameters defining resolution and aspect ratio (calculated from bounds)
   */ 
  discretized_line_segment(const Eigen::Ref<const Vector2> &zeta_s, // start point in linear coords
                           const Eigen::Ref<const Vector2> &zeta_e, // end point in linear coords
                           std::shared_ptr<physical_parameters> phys,
                           const simulation_parameters &sim);

  /*! evaluate linear coords at point along segment
   * \param lam parameter along segment in [0,1]
   * \return point at lam along segment in linear coords
   */
  Vector2 get_z(const T &lam) const;

  /*! evaluate s-p-coords at point along segment
   * \param lam parameter along segment in [0,1]
   * \return point at lam along segment in s-p-coords
   */
  Vector2 get_x(const T &lam) const;

  /*! evaluate s-p-coord tangent at point along segment
   * \param lam parameter along segment in [0,1]
   * \return tangent at point at lam along segment in s-p-coords
   */
  Vector2 get_t(const T &lam) const;

  /*! increase / set resolution of segment
   * \param newlen new number of points to use, if negative double number of segments clamped by max_resolution
   */
  void refine(int newlen = -1);

  /*! return number of points used for discretization */
  int size() const;
  
  /*! return incomplete string representation (mainly for python) */
  std::string repr() const;

  /*! bind to python module */
  static void bind(py::module_ &m);
};

#endif
