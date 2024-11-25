
#include "discretized_line_segment2.hpp"

template class discretized_line_segment<double>;
template class discretized_line_segment<long double>;

#define DLS discretized_line_segment<T>

template <typename T>
DLS::discretized_line_segment(const Eigen::Ref<const DLS::Vector2> &zeta_s,
                              const Eigen::Ref<const DLS::Vector2> &zeta_e,
                              const physical_parameters &phys,
                              const simulation_parameters &sim) :
  start(zeta_s), end(zeta_e), direction((zeta_e - zeta_s).normalized()),
  phys(phys), sim(sim),
  z(sim.min_line_resolution, 2),
  x(sim.min_line_resolution, 2),
  t(sim.min_line_resolution, 2),
  s(sim.min_line_resolution - 1, 2),
  errb(0), area(0) {
  
  for (int i = 0; i < sim.min_line_resolution; ++i) {
    T lam = i / (T)(sim.min_line_resolution - 1);
    z.row(i) = get_z(lam);
    x.row(i) = get_x(lam);
    t.row(i) = get_t(lam);
  }
  for (int i = 1; i < sim.min_line_resolution; ++i) {
    // calculate intersection point
    // TODO: find a more numerically stable (and maybe faster) way to calculate this [also in refine]
    T c = t.row(i).cross(t.row(i-1));
    T l = (c == 0) ? (T)0.5 : std::clamp(t.row(i).cross(x.row(i) - x.row(i-1)) / c, (T)0, (T)1);
    s.row(i-1) = x.row(i-1) + l * t.row(i-1);
    
    // calculate area
    area += 0.25 * ((x(i,0) - x(i-1,0)) * (x(i,1) + x(i-1,1))
                    + (s(i-1,0) - x(i-1,0)) * (s(i-1,1) + x(i-1,1))
                    + (x(i,0) - s(i-1,0)) * (x(i,1) + s(i-1,1)));
    
    // calculate error bound
    errb += 0.25 * std::abs((x.row(i) - x.row(i-1)).cross(s.row(i-1) - x.row(i-1)));
  }

  while (size() < sim.max_line_resolution && errb > sim.line_tolerance)
    refine();
}

template <typename T>
DLS::Vector2 DLS::get_z(const T &lam) const {
  return start * (1 - lam) + end * lam;
}

template <typename T>
DLS::Vector2 DLS::get_x(const T &lam) const {
  return phys.itf<T>(get_z(lam));
}

template <typename T>
DLS::Vector2 DLS::get_t(const T &lam) const {
  return phys.ditf<T>(get_z(lam)).array() * direction.array();
}

template <typename T>
void DLS::refine() {
  // new length: adding a point in the middle of every gap
  int newlen = size() * 2 - 1;

  // make new arrays, reset values
  Matrix2X newz(newlen, 2);
  Matrix2X newx(newlen, 2);
  Matrix2X newt(newlen, 2);
  s = Matrix2X(newlen - 1, 2);
  errb = 0;
  area = 0;

  // create new / copy old points
  for (int i = 0; i < newlen; ++i) {
    T lam = i / (T)(newlen - 1);

    // calculate new points and copy old ones
    if ((i & 1)) {
      newz.row(i) = get_z(lam);
      newx.row(i) = get_x(lam);
      newt.row(i) = get_t(lam);
    } else {
      newz.row(i) = z.row(i / 2);
      newx.row(i) = x.row(i / 2);
      newt.row(i) = t.row(i / 2);
    }
  }
  // recalc area and errb
  for (int i = 1; i < newlen; ++i) {
      // calculate intersection point
      // TODO: find a more numerically stable way to calculate this
      T c = newt.row(i).cross(newt.row(i-1));
      T l = (c == 0) ? (T)0.5 : std::clamp(newt.row(i).cross(newx.row(i) - newx.row(i-1)) / c, (T)0, (T)1);
      s.row(i-1) = newx.row(i-1) + l * newt.row(i-1);

      // calculate area
      area += 0.25 * ((newx(i,0) - newx(i-1,0)) * (newx(i,1) + newx(i-1,1))
                      + (s(i-1,0) - newx(i-1,0)) * (s(i-1,1) + newx(i-1,1))
                      + (newx(i,0) - s(i-1,0)) * (newx(i,1) + s(i-1,1)));
      
      // calculate error bound
      errb += 0.25 * std::abs((newx.row(i) - newx.row(i-1)).cross(s.row(i-1) - newx.row(i-1)));
  }

  z = std::move(newz);
  x = std::move(newx);
  t = std::move(newt);
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
    .def_readonly("z", &DLS::z)
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
