#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


void hello() {
  py::print("Hello from C++!");
}

PYBIND11_MODULE(_atmosphere_bgs, m) {
  m.def("_hello", &hello);
}
