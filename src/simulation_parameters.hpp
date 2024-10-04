#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include "common.hpp"

struct simulation_parameters
{
  // limits for s and p
  Eigen::Vector2d spmin, spmax;

  int boundary_res;

  double area_tolerance;
  int max_refine_steps;

  double line_tolerance;
  int min_line_resolution;
  int max_line_resolution;
  
  simulation_parameters(double smin = 0.1747087883522576,
                        double smax = 0.9809294113733709,
                        double pmin = 7300,
                        double pmax = 141855,
                        int boundary_res = 1000,
                        double area_tolerance = 1e-4,
                        int max_refine_steps = 1000,
                        double line_tolerance = 1e-3,
                        int min_line_resolution = 4,
                        int max_line_resolution = 100) :
    spmin(smin, pmin), spmax(smax, pmax), boundary_res(boundary_res),
    area_tolerance(area_tolerance), max_refine_steps(max_refine_steps),
    line_tolerance(line_tolerance), min_line_resolution(min_line_resolution), max_line_resolution(max_line_resolution) {}
  
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
       << ", max_line_resolution = " << max_line_resolution;
    return ss.str();
  }

  static void bind(py::module_ &m) {
    py::class_<simulation_parameters>(m, "SimulationParameters")
      .def(py::init<double, double, double, double,
           int, double, int, double, int, int>(),
           py::arg("smin") = 0.1747087883522576, py::arg("smax") = 0.9809294113733709,
           py::arg("pmin") = 7300., py::arg("pmax") = 141855.,
           py::arg("boundary_res") = 1000,
           py::arg("area_tolerance") = 1e-4,
           py::arg("max_refine_steps") = 1000,
           py::arg("line_tolerance") = 1e-3,
           py::arg("min_line_resolution") = 4,
           py::arg("max_line_resolution") = 100)
      .def_readonly("spmin", &simulation_parameters::spmin)
      .def_readonly("spmax", &simulation_parameters::spmax)
      .def_readonly("boundary_res", &simulation_parameters::boundary_res)
      .def_readonly("area_tolerance", &simulation_parameters::area_tolerance)
      .def_readonly("max_refine_steps", &simulation_parameters::max_refine_steps)
      .def_readonly("line_tolerance", &simulation_parameters::line_tolerance)
      .def_readonly("min_line_resolution", &simulation_parameters::min_line_resolution)
      .def_readonly("max_line_resolution", &simulation_parameters::max_line_resolution)
      .def("__repr__", &simulation_parameters::repr);
  }
};

#endif
