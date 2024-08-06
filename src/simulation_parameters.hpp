#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include "common.hpp"

struct simulation_parameters
{
  // limits for s and p
  Eigen::Vector2d spmin, spmax;

  int boundary_res;

  double tol;
  int max_refine_steps;
  
  simulation_parameters(double smin = 0.1747087883522576,
                        double smax = 0.9809294113733709,
                        double pmin = 7300,
                        double pmax = 141855,
                        int boundary_res = 1000,
                        double tol = 1e-4,
                        int max_refine_steps = 100) :
    spmin(smin, pmin), spmax(smax, pmax), boundary_res(boundary_res), tol(tol), max_refine_steps(max_refine_steps) {}
  
  std::string repr() const {
    std::stringstream ss;
    ss << "smin = " << spmin(0)
       << ", smax = " << spmax(0)
       << ", pmin = " << spmin(1)
       << ", pmax = " << spmax(1)
       << ", boundary_res = " << boundary_res
       << ", tol = " << tol
       << ", max_refine_steps = " << max_refine_steps;
    return ss.str();
  }
};

#define BIND_SIMULATION_PARAMETERS(m)                                   \
  py::class_<simulation_parameters>(m, "SimulationParameters")          \
  .def(py::init<double, double, double, double, int, double, int>(),    \
       py::arg("smin") = 0.1747087883522576, py::arg("smax") = 0.9809294113733709, \
       py::arg("pmin") = 7300., py::arg("pmax") = 141855.,              \
       py::arg("boundary_res") = 1000, py::arg("tol") = 0.1,            \
       py::arg("max_refine_steps") = 100)                               \
  .def_readonly("spmin", &simulation_parameters::spmin)                 \
  .def_readonly("spmax", &simulation_parameters::spmax)                 \
  .def_readonly("boundary_res", &simulation_parameters::boundary_res)   \
  .def_readonly("tol", &simulation_parameters::tol)                     \
  .def_readonly("max_refine_steps", &simulation_parameters::max_refine_steps) \
  .def("__repr__", &simulation_parameters::repr)

#endif
