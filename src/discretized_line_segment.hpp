#ifndef DISCRETIZED_LINE_SEGMENT
#define DISCRETIZED_LINE_SEGMENT

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>


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
};

#define BIND_TANGENT_POINT(m)                  \
  py::class_<tangent_point>(m, "TangentPoint") \
  .def("__repr__", &tangent_point::repr)       \
  .def_readonly("zeta", &tangent_point::zeta)  \
  .def_readonly("x", &tangent_point::x)        \
  .def_readonly("t", &tangent_point::t)        \
  .def_readonly("s", &tangent_point::s);

struct discretized_line_segment
{
private:
  double sqr(const double &x) const { return x * x; }
  
public:
  physical_parameters phys;
  std::list<tangent_point> points;
  Eigen::Vector2d dzeta;
  
  std::multimap<double, typename std::list<tangent_point>::iterator> errmap;
  double errb;
  double area;
    
  discretized_line_segment(const Eigen::Ref<const Eigen::Vector2d> &zeta_s, // start point in linear coords
                           const Eigen::Ref<const Eigen::Vector2d> &zeta_e, // end point in linear coords
                           const physical_parameters &phys,
                           const double &tol) :
    phys(phys),
    dzeta(zeta_e - zeta_s) {

    const int min_res = 2;
    const int max_res = 100;
    
    tangent_point ts(zeta_s, dzeta, this);
    tangent_point te(zeta_e, dzeta, this);
    ts.update_intersect(te);
      
    points.push_back(ts);
    points.push_back(te);

    errb = points.front().errb(points.back());
    errmap.insert({-errb, points.begin()});

    area = ts.area(te);
    
    for (int i = 0; (i < max_res) && (i < min_res || errb > tol); ++i)
      refine();
  }

  void refine() {
    auto f = *errmap.begin();
    errmap.erase(errmap.begin());
    
    // std::cout << "errb = " << errb;
    
    auto it0 = f.second++;
    auto it2 = f.second;

    // todo: more efficient way to calculate area
    area -= it0->area(*it2);
    
    auto it1 = points.insert(it2, tangent_point(0.5 * (it0->zeta + it2->zeta), dzeta, this));
    
    it0->update_intersect(*it1);
    it1->update_intersect(*it2);
    
    double err0 = it0->errb(*it1);
    double err1 = it1->errb(*it2);
    
    errb += f.first + err0 + err1;

    // todo: more efficient way to calculate area (just add new triangle)
    area += it0->area(*it1) + it1->area(*it2);
    
    // std::cout << " replacing with " << err0  << " + " << err1 << std::endl;
    
    errmap.insert({-err0, it0});
    errmap.insert({-err1, it1});
  }

  // get (a) tangent at point x on the line
  inline Eigen::Vector2d get_tangent(const Eigen::Ref<const Eigen::Vector2d> &zeta,
                              const Eigen::Ref<const Eigen::Vector2d> &dzeta) const {
    return phys.ditf(zeta).array() * dzeta.array();
  }

  double length() const {
    double res = 0;
    for (auto it = points.begin(), itl = it++; it != points.end(); itl = it++) {
      res += (it->x - itl->x).norm();
    }
    return res;
  }
  
  std::string repr() const {
    std::stringstream ss;
    ss << "DiscretizedLineSegment(length = " << points.size() << ")";
    return ss.str();
  }
};

#define BIND_DISCRETIZED_LINE_SEGMENT(m)                             \
  py::class_<discretized_line_segment>(m, "DiscretizedLineSegment")  \
  .def(py::init<const Eigen::Ref<const Eigen::Vector2d> &,           \
       const Eigen::Ref<const Eigen::Vector2d> &,                    \
       const physical_parameters &, const double &>())               \
  .def("__repr__", &discretized_line_segment::repr)                  \
  .def_readonly("area", &discretized_line_segment::area)             \
  .def("refine", &discretized_line_segment::refine)                  \
  .def_readonly("errb", &discretized_line_segment::errb)             \
  .def_readonly("points", &discretized_line_segment::points)         \
  .def("length", &discretized_line_segment::length);


void tangent_point::update_intersect(const tangent_point &next) {
  double c = next.t.cross(t), lam;
  if (c == 0) {
    lam = 0.5;
  } else {
    lam = next.t.cross(next.x - x) / c;
    // double lam = (next.t0 * (next.x1 - x1) - next.t1 * (next.x0 - x0)) / (next.t0 * t1 - next.t1 * t0);
  }
  s = x + lam * t;
}

double tangent_point::errb(const tangent_point &next) const {
  return 0.25 * std::abs((next.x - x).cross(s - x));
}

double tangent_point::area(const tangent_point &next) const {
  return 0.25 * ((next.x(0) - x(0)) * (next.x(1) + x(1))
                 + (s(0) - x(0)) * (s(1) + x(1))
                 + (next.x(0) - s(0)) * (next.x(1) + s(1)));
}

tangent_point::tangent_point(const Eigen::Ref<const Eigen::Vector2d> &zeta,
                             const Eigen::Ref<const Eigen::Vector2d> &dzeta,
                             const discretized_line_segment *ls) :
  zeta(zeta), x(ls->phys.itf(zeta)), t(ls->get_tangent(zeta, dzeta)),
  s(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()) {
}

std::string tangent_point::repr() const {
  std::stringstream ss;
  ss << "TangentPoint("
     << zeta(0) << ", " << zeta(1) << "; "
     << x(0) << ", " << x(1) << "; "
     << t(0) << ", " << t(1) << "; "
     << s(0) << ", " << s(1) << ")";
  return ss.str();
}

#endif
