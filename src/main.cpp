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
  BIND_TANGENT_POINT(m);
  BIND_DISCRETIZED_LINE_SEGMENT(m); 

  BIND_PDEDGE(m);
  BIND_PDMESH(m);
  BIND_HALFSPACE_INTERSECTION(m);

  BIND_DIAGRAM_EDGE(m);
  BIND_LAGUERRE_DIAGRAM(m);

  BIND_RASTERIZER_SEGMENT(m);
  BIND_RASTERIZER_EVENT(m);
  BIND_RASTERIZER(m);
}
