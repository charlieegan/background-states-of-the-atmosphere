
#include "discretized_line_segment.hpp"

template class discretized_line_segment<double>;
template class discretized_line_segment<long double>;

#define DLS discretized_line_segment<T>

template <typename T>
DLS::discretized_line_segment(const Eigen::Ref<const DLS::Vector2> &zeta_s,
                              const Eigen::Ref<const DLS::Vector2> &zeta_e,
                              const physical_parameters &phys,
                              const simulation_parameters &sim) :
  start(zeta_s), end(zeta_e), direction((zeta_e - zeta_s).normalized()),
  phys(&phys), 
  max_resolution(sim.max_line_resolution), aspect(sim.spmax(1) - sim.spmin(1), sim.spmax(0) - sim.spmin(0)),
  lams(sim.min_line_resolution), x(sim.min_line_resolution, 2),
  errb(0), area(0) {

  // handle length 0 
  if (zeta_s == zeta_e) {
    lams.resize(2);
    x.resize(2, 2);
    t.resize(2, 2);
    s.resize(2 - 1, 2);

    lams(0) = 0;
    lams(1) = 1;
    s.row(0) = x.row(0) = x.row(1) = get_x(0);
    t.array() = 0;
    return;
  }
  
  // calculate initial positions
  for (int i = 0; i < sim.min_line_resolution; ++i) {
    T lam = i / (T)(sim.min_line_resolution - 1);
    lams(i) = lam;
    x.row(i) = get_x(lam);
  }

  // refine (this also calculates t and s)
  refine(sim.min_line_resolution);

  // continue refining while tolerance is not reached
  while (size() < sim.max_line_resolution && errb > sim.line_tolerance)
    refine();
}

template <typename T>
DLS::Vector2 DLS::get_z(const T &lam) const {
  return start * (1 - lam) + end * lam;
}

template <typename T>
DLS::Vector2 DLS::get_x(const T &lam) const {
  return phys->itf<T>(get_z(lam));
}

template <typename T>
DLS::Vector2 DLS::get_t(const T &lam) const {
  return phys->ditf<T>(get_z(lam)).array() * direction.array();
}

template <typename T>
void DLS::refine(int newlen) {
  // don't refine length 0 segment (since it is already exact and would lead to NaNs)
  if (start == end)
    return;
  
  int oldlen = size();
  
  // if no new length is given, double current (split every current segment in half)
  if (newlen < 2) {
    newlen = std::min(oldlen * 2 - 1, max_resolution);
    if (newlen <= oldlen)
      return;
  }

  // calculate path length
  VectorX lens(oldlen);
  lens(0) = 0;
  for (int i = 1; i < oldlen; ++i) {
    T len = ((x.row(i) - x.row(i - 1)).transpose().array() * aspect.array()).matrix().norm();
    lens(i) = lens(i - 1) + len;
  }
  T totlen = lens(oldlen - 1);

  // calculate new lam positions to be approx evenly spaced in s-p-coords
  VectorX oldlams = lams;
  lams.resize(newlen);
  for (int i = 0, j = 1; j < newlen - 1; ++j) {
    T oj = totlen * j / (newlen - 1); // evenly spaced position along length
    while (i + 2 < oldlen && oj > lens(i + 1)) ++i; // find index i s.t. lens[i] <= oj < lens[i+1]
    T sl = (oldlams(i + 1) - oldlams(i)) / (lens(i + 1) - lens(i)); // (inverse) slope of length-line-segment
    lams(j) = (oj - lens(i)) * sl + oldlams(i);
  }
  lams(0) = (T)0.0;
  lams(newlen - 1) = (T)1.0;
  
  // make new resized arrays, reset values
  x.resize(newlen, 2);
  t.resize(newlen, 2);
  s.resize(newlen - 1, 2);
  errb = 0;
  area = 0;

  // calculate x, t for new points
  for (int i = 0; i < newlen; ++i) {
    T lam = lams(i);
    x.row(i) = get_x(lam);
    t.row(i) = get_t(lam);
  }
  // calculate s, area and errb for new points
  for (int i = 1; i < newlen; ++i) {
      // calculate intersection point
      // TODO: find a more numerically stable way to calculate this
      T c = t.row(i).cross(t.row(i-1));
      T u = std::abs(t.row(i).dot(x.row(i) - x.row(i-1)) / t.row(i).dot(t.row(i)));
      T l = (c == 0) ? (T)0.5 * u : std::clamp(t.row(i).cross(x.row(i) - x.row(i-1)) / c, (T)0, u);
      s.row(i-1) = x.row(i-1) + l * t.row(i-1);

      // auto b = x.row(i) - x.row(i-1);
      // T tan1 = t.row(i-1).cross(b) / t.row(i-1).dot(b);
      // T tan2 = -t.row(i).cross(b) / t.row(i).dot(b);
      // T e2 = 0.5 * b.dot(b) * std::abs(tan1 * tan2 / (tan1 + tan2));

      // area of tangent / point triangle
      T e = 0.5 * std::abs((x.row(i) - x.row(i-1)).cross(s.row(i-1) - x.row(i-1)));

      // if (std::abs(std::copysign(e2, e) - e) > 1e-9)
      //   py::print(FORMAT("e = {}, e2 = {}, tan1 = {}, tan2 = {}, x^2 = {}", e, e2, tan1, tan2, b.dot(b)));
      
      // calculate area
      area += 0.5 * ((x(i,0) - x(i-1,0)) * (x(i,1) + x(i-1,1))) + 0.5 * e * std::copysign(e, x(i,0) - x(i-1,0));
      
      // area += 0.25 * ((x(i,0) - x(i-1,0)) * (x(i,1) + x(i-1,1))
      //                 + (s(i-1,0) - x(i-1,0)) * (s(i-1,1) + x(i-1,1))
      //                 + (x(i,0) - s(i-1,0)) * (x(i,1) + s(i-1,1)));
      
      // calculate error bound
      errb += 0.5 * e;
  }

  //py::print(FORMAT("res {}, area {}, errb {}", newlen, area, errb));
}

template <typename T>
int DLS::size() const {
  return x.rows();
}

template <typename T>
std::string DLS::repr() const {
  return FORMAT("DiscretizedLineSegment(size = {})", size());
}

template <typename T>
void DLS::bind(py::module_ &m) {
  py::class_<DLS>(m, ("DiscretizedLineSegment_" + type_name<T>::value()).c_str())
    .def(py::init<const Eigen::Ref<const DLS::Vector2> &,
         const Eigen::Ref<const DLS::Vector2> &,
         const physical_parameters &,
         const simulation_parameters &>())
    .def_readonly("start", &DLS::start)
    .def_readonly("end", &DLS::end)
    .def_readonly("direction", &DLS::direction)
    // .def_readonly("z", &DLS::z)
    .def_readonly("lams", &DLS::lams)
    .def_readonly("x", &DLS::x)
    .def_readonly("t", &DLS::t)
    .def_readonly("s", &DLS::s)
    .def("get_z", &DLS::get_z)
    .def("get_x", &DLS::get_x)
    .def("get_t", &DLS::get_t)
    .def("refine", &DLS::refine)
    .def("size", &DLS::size)
    .def("__repr__", &DLS::repr)
    .def_readonly("area", &DLS::area)
    .def_readonly("errb", &DLS::errb);

}

#undef DLS
