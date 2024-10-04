#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#define PROFILING

#define DEBUG_CHECKS
#include "rasterizer.hpp"
#undef DEBUG_CHECKS

#include "laguerre_diagram.hpp"
#include "timer.hpp"
#include "halfspace_intersection.hpp"


void hello() {
  py::print("Hello from C++!"); 
}
 

PYBIND11_MODULE(_atmosphere_bgs, m) {
  m.def("_hello", &hello); 

  BIND_TIMER(m);
  
  BIND_PHYSICAL_PARAMETERS(m);
  BIND_SIMULATION_PARAMETERS(m);

  tangent_point::bind(m);
  discretized_line_segment::bind(m);

  BIND_PDEDGE(m);
  BIND_PDMESH(m, double, "PDMesh_f64");
  BIND_HALFSPACE_INTERSECTION(m, double, "HalfspaceIntersection_fp64");

  BIND_PDMESH(m, long double, "PDMesh_f128");
  BIND_HALFSPACE_INTERSECTION(m, long double, "HalfspaceIntersection_fp128");

  laguerre_diagram<double>::diagram_edge::bind(m);
  laguerre_diagram<double>::bind(m);
  laguerre_diagram<long double>::diagram_edge::bind(m);
  laguerre_diagram<long double>::bind(m);

  rasterizer::segment::bind(m);
  rasterizer::event::bind(m);
  rasterizer::bind(m);
}
