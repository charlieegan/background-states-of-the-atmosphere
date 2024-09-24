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
  BIND_PDMESH(m, double, "PDMesh_f64");
  BIND_HALFSPACE_INTERSECTION(m, double, "HalfspaceIntersection_fp64");

  BIND_DIAGRAM_EDGE(m, double, "DiagramEdge_f64");
  BIND_LAGUERRE_DIAGRAM(m, double, "LaguerreDiagram_f64");


  BIND_PDMESH(m, long double, "PDMesh_f128");
  BIND_HALFSPACE_INTERSECTION(m, long double, "HalfspaceIntersection_fp128");

  BIND_DIAGRAM_EDGE(m, long double, "DiagramEdge_f128");
  BIND_LAGUERRE_DIAGRAM(m, long double, "LaguerreDiagram_f128");

  m.def("LaguerreDiagram",
        [](const laguerre_diagram<double>::seeds_t &ys,
           const Eigen::Ref<const laguerre_diagram<double>::VecX> &duals,
           const physical_parameters &phys,
           const simulation_parameters &sim) {
          return laguerre_diagram<double>(ys, duals, phys, sim);
        },
        py::arg("ys"), py::arg("duals"),          
        py::arg("phys"), py::arg("sim"));

  m.def("LaguerreDiagram",
        [](const laguerre_diagram<long double>::seeds_t &ys,
           const Eigen::Ref<const laguerre_diagram<long double>::VecX> &duals,
           const physical_parameters &phys,
           const simulation_parameters &sim) {
          return laguerre_diagram<long double>(ys, duals, phys, sim);
        },
        py::arg("ys"), py::arg("duals"),          
        py::arg("phys"), py::arg("sim"));

  BIND_RASTERIZER_SEGMENT(m);
  BIND_RASTERIZER_EVENT(m);
  BIND_RASTERIZER(m);
}
