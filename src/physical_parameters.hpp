#ifndef PHYSICAL_PARAMETERS_HPP
#define PHYSICAL_PARAMETERS_HPP

#include "common.hpp"

struct physical_parameters
{
  double a, Omega, p00, kappa, cp;
  double ikappa, tc0p, tc1p;

  physical_parameters(double a = 6371000.,
                      double Omega = 7.2921e-05,
                      double p00 = 101325.,
                      double kappa = 2. / 7.,
                      double cp = 1003.5) :
    a(a), Omega(Omega), p00(p00), kappa(kappa), cp(cp),
    ikappa(1. / kappa), tc0p(1. / (a * a)), tc1p(kappa * cp * std::pow(p00, -kappa)) {}

  template <typename T>
  Eigen::Vector2<T> itf(const Eigen::Ref<const Eigen::Vector2<T>> &zeta) const {
    return Eigen::Vector2<T>(std::sqrt(std::max((T)0., (T)1. - 1. / zeta(0))),
                             p00 * std::pow(std::max((T)0., zeta(1)), ikappa));
  }

  template <typename T>
  Eigen::Vector2<T> ditf(const Eigen::Ref<const Eigen::Vector2<T>> &zeta) const {
    return Eigen::Vector2<T>(0.5 / (std::sqrt(std::max((T)0., (T)1. - 1. / zeta(0))) * zeta(0) * zeta(0)),
                             p00 * ikappa * std::pow(std::max((T)0., (T)zeta(1)), ikappa - 1));
  }

  template <typename T>
  Eigen::Vector2<T> tf(const Eigen::Ref<const Eigen::Vector2<T>> &x) const {
    return Eigen::Vector2<T>(1. / (1. - x[0] * x[0]),
                             std::pow(x[1] / p00, kappa));
  }

  template <typename T>
  Eigen::Vector2<T> dtf(const Eigen::Ref<const Eigen::Vector2<T>> &x) const {
    return Eigen::Vector2<T>(2. * x[0] / ((1. - x[0] * x[0]) * (1. - x[0] * x[0])),
                             kappa / std::pow(p00, kappa) * std::pow(x[1], kappa - 1));
  }

  std::string repr() const {
    std::stringstream ss;
    ss << "a = " << a
       << ", Omega = " << Omega
       << ", p00 = " << p00
       << ", kappa = " << kappa
       << ", cp = " << cp;
    return ss.str();
  }

  static void bind(py::module_ &m) {
    py::class_<physical_parameters>(m, "PhysicalParameters")
      .def(py::init<double, double, double, double, double>(),
           py::arg("a") = 6371000., py::arg("Omega") = 7.2921e-05,
           py::arg("p00") = 101325., py::arg("kappa") = 2. / 7., py::arg("cp") = 1003.5)
      .def("__repr__", &physical_parameters::repr)
      .def("itf", &physical_parameters::itf<double>, py::arg("zeta"))
      .def("tf", &physical_parameters::tf<double>, py::arg("x"))
      .def_readonly("a", &physical_parameters::a)
      .def_readonly("Omega", &physical_parameters::Omega)
      .def_readonly("p00", &physical_parameters::p00)
      .def_readonly("kappa", &physical_parameters::kappa)
      .def_readonly("cp", &physical_parameters::cp);
  }
};

#endif
