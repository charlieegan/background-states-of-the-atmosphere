/*! \mainpage Background states of the atmosphere
 *
 * This page contains documentation for the code generated with doxygen.
 * Currently the documentation is very bare-bones but this will hopefully improve in the future.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

namespace py = pybind11;

#define PROFILING

//#define EXPENSIVE_DEBUG_CHECKS
//#define DEBUG_CHECKS
#include "rasterizer.hpp"
//#undef DEBUG_CHECKS

#include "laguerre_diagram.hpp"
#include "timer.hpp"
#include "halfspace_intersection.hpp"
 
/// @private
PYBIND11_MODULE(_atmosphere_bgs, m) {
  timer::bind(m);

  physical_parameters::bind(m);
  simulation_parameters::bind(m);
  
  discretized_line_segment<double>::bind(m);

  pdedge::bind(m);
  pdmesh<double>::bind(m);
  pdmesh<long double>::bind(m);
  halfspace_intersection<double>::bind(m);
  halfspace_intersection<long double>::bind(m);

  laguerre_diagram<double>::diagram_edge::bind(m);
  laguerre_diagram<double>::bind(m);
  laguerre_diagram<long double>::diagram_edge::bind(m);
  laguerre_diagram<long double>::bind(m);

  rasterizer::segment::bind(m);
  rasterizer::event::bind(m);
  rasterizer::bind(m);
}
