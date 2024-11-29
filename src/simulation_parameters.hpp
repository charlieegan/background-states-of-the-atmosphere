#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include "common.hpp"

struct simulation_parameters
{
  // limits for s and p
  Eigen::Vector2d spmin, spmax;

  // resolution of top cell boundary
  int boundary_res;

  // relative error to try to reach with discretization
  double area_tolerance;
  
  // maximal number of times to try to reach error (does not do much)
  int max_refine_steps;

  // relative area error to try to reach per line
  double line_tolerance;

  // minimal number of points to use per line 
  int min_line_resolution;

  // maximal number of points to use per line
  int max_line_resolution;

  // factor to use for negative areas (set <= 0 to not allow negative areas)
  double negative_area_scaling;
  
  simulation_parameters(double smin = 0.1747087883522576,
                        double smax = 0.9809294113733709,
                        double pmin = 7300,
                        double pmax = 141855,
                        int boundary_res = 1000,
                        double area_tolerance = 1e-4,
                        int max_refine_steps = 1000,
                        double line_tolerance = 1e-3,
                        int min_line_resolution = 4,
                        int max_line_resolution = 100,
                        double negative_area_scaling = 0.0) :
    spmin(smin, pmin), spmax(smax, pmax), boundary_res(boundary_res),
    area_tolerance(area_tolerance), max_refine_steps(max_refine_steps),
    line_tolerance(line_tolerance), min_line_resolution(min_line_resolution), max_line_resolution(max_line_resolution),
    negative_area_scaling(negative_area_scaling) {}

  void set_smin(const double &v) { spmin[0] = v; }
  void set_smax(const double &v) { spmax[0] = v; }
  void set_pmin(const double &v) { spmin[1] = v; }
  void set_pmax(const double &v) { spmax[1] = v; }
  double get_smin() const { return spmin[0]; }
  double get_smax() const { return spmax[0]; }
  double get_pmin() const { return spmin[1]; }
  double get_pmax() const { return spmax[1]; }
  
  std::string repr() const {
    std::stringstream ss;
    ss << "smin = " << spmin(0)
       << ", smax = " << spmax(0)
       << ", pmin = " << spmin(1)
       << ", pmax = " << spmax(1)
       << ", boundary_res = " << boundary_res
       << ", area_tolerance = " << area_tolerance
       << ", max_refine_steps = " << max_refine_steps
       << ", line_tolerance = " << line_tolerance
       << ", min_line_resolution = " << min_line_resolution
       << ", max_line_resolution = " << max_line_resolution
       << ", negative_area_scaling = " << negative_area_scaling;
    return ss.str();
  }

  static py::tuple pickle(const simulation_parameters &s) {
    return py::make_tuple(s.spmin[0], s.spmax[0], s.spmin[1], s.spmax[1], 
                          s.boundary_res, s.area_tolerance, s.max_refine_steps,
                          s.line_tolerance, s.min_line_resolution, s.max_line_resolution,
                          s.negative_area_scaling);
  }

  static simulation_parameters unpickle(py::tuple t) {
    if (t.size() != 11)
      throw std::runtime_error("invalid state in simulation_parameters::unpickle");

    return simulation_parameters(t[0].cast<double>(),
                                 t[1].cast<double>(),
                                 t[2].cast<double>(),
                                 t[3].cast<double>(),
                                 t[4].cast<int>(),
                                 t[5].cast<double>(), 
                                 t[6].cast<int>(), 
                                 t[7].cast<double>(), 
                                 t[8].cast<int>(), 
                                 t[9].cast<int>(),
                                 t[10].cast<double>());
  }
  
  static void bind(py::module_ &m) {
    py::class_<simulation_parameters>(m, "SimulationParameters")
      .def(py::init<double, double, double, double,
           int, double, int, double, int, int, double>(),
           py::arg("smin") = 0.1747087883522576, py::arg("smax") = 0.9809294113733709,
           py::arg("pmin") = 7300., py::arg("pmax") = 141855.,
           py::arg("boundary_res") = 1000,
           py::arg("area_tolerance") = 1e-4,
           py::arg("max_refine_steps") = 1000,
           py::arg("line_tolerance") = 1e-3,
           py::arg("min_line_resolution") = 4,
           py::arg("max_line_resolution") = 100,
           py::arg("negative_area_scaling") = 0.0)
      .def_readwrite("spmin", &simulation_parameters::spmin)
      .def_readwrite("spmax", &simulation_parameters::spmax)
      .def_property("smin", &simulation_parameters::get_smin, &simulation_parameters::set_smin)
      .def_property("smax", &simulation_parameters::get_smax, &simulation_parameters::set_smax)
      .def_property("pmin", &simulation_parameters::get_pmin, &simulation_parameters::set_pmin)
      .def_property("pmax", &simulation_parameters::get_pmax, &simulation_parameters::set_pmax)
      .def_readwrite("boundary_res", &simulation_parameters::boundary_res)
      .def_readwrite("area_tolerance", &simulation_parameters::area_tolerance)
      .def_readwrite("max_refine_steps", &simulation_parameters::max_refine_steps)
      .def_readwrite("line_tolerance", &simulation_parameters::line_tolerance)
      .def_readwrite("min_line_resolution", &simulation_parameters::min_line_resolution)
      .def_readwrite("max_line_resolution", &simulation_parameters::max_line_resolution)
      .def_readwrite("negative_area_scaling", &simulation_parameters::negative_area_scaling)
      .def(py::pickle(&simulation_parameters::pickle, &simulation_parameters::unpickle))
      .def("__repr__", &simulation_parameters::repr);
  }
};

#endif
