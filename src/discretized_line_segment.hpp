#ifndef DISCRETIZED_LINE_SEGMENT
#define DISCRETIZED_LINE_SEGMENT

#include "common.hpp"

struct discretized_line_segment;

// a class representing a point on an inner boundary of a laguerre cell
struct tangent_point
{
  // linearized coords (zeta)
  // original (s-p-) coords (x)
  // tangent in s-p-coords (t)
  // intersection of this tangent with next tangent (s)
  Eigen::Vector2d zeta, x, t, s;
  
  tangent_point() {}
  tangent_point(const tangent_point &o) :
    zeta(o.zeta), x(o.x), t(o.t), s(o.s) {}
      
  // construct given linearized coords zeta0, zeta1 and
  // intersect will not be calculated in this constructor
  tangent_point(const Eigen::Ref<const Eigen::Vector2d> &zeta,
                const Eigen::Ref<const Eigen::Vector2d> &dzeta,
                const discretized_line_segment *ls);

  // update the tangent intersection point of this with next
  void update_intersect(const tangent_point &next);

  // find the error bound for segment between this and next (requires this->is* to be up to date)
  double errb(const tangent_point &next) const;

  double area(const tangent_point &next) const;
    
  std::string repr() const;

  static void bind(py::module_ &m);
};

struct discretized_line_segment
{
public:
  const physical_parameters &phys;
  std::list<tangent_point> points;
  Eigen::Vector2d dzeta;
  
  std::multimap<double, typename std::list<tangent_point>::iterator> errmap;
  double errb;
  double area;
    
  discretized_line_segment(const Eigen::Ref<const Eigen::Vector2d> &zeta_s, // start point in linear coords
                           const Eigen::Ref<const Eigen::Vector2d> &zeta_e, // end point in linear coords
                           const physical_parameters &phys,
                           const simulation_parameters &sim);

  void refine();

  // get (a) tangent at point x on the line
  Eigen::Vector2d get_tangent(const Eigen::Ref<const Eigen::Vector2d> &zeta,
                                     const Eigen::Ref<const Eigen::Vector2d> &dzeta) const;

  double length() const;
  
  std::string repr() const;

  static void bind(py::module_ &m);
};

#endif
