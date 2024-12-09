#ifndef PHYSICAL_PARAMETERS_HPP
#define PHYSICAL_PARAMETERS_HPP

#include "common.hpp"

struct physical_parameters
{
  double a;      //!< earth radius [m]
  double Omega;  //!< earth angular velocity [rad/s]
  double p00;    //!< standard atmosphere pressure [Pa]
  double kappa;  //!< 
  double cp;     //!< specific heat of air [J/kgK] ?
  
  double ikappa; //!< 1 / kappa (for more efficient internal usage)

  physical_parameters(double a = 6371000.,
                      double Omega = 7.2921e-05,
                      double p00 = 101325.,
                      double kappa = 2. / 7.,
                      double cp = 1003.5) :
    a(a), Omega(Omega), p00(p00), kappa(kappa), cp(cp),
    ikappa(1. / kappa) {}

  /*!
   * Transform a point linearized coords -> physical s-p-coords
   * \param zeta linear coord point
   * \return zeta transformed to s-p-coords
   */
  template <typename T>
  Eigen::Vector2<T> itf(const Eigen::Ref<const Eigen::Vector2<T>> &zeta) const {
    return Eigen::Vector2<T>(std::sqrt(std::max((T)0., (T)1. - 1. / zeta(0))),
                             p00 * std::pow(std::max((T)0., zeta(1)), ikappa));
  }

  /*!
   * Derivative (diagonal of Jacobian, off-diagonal is 0) of transform linearized coords -> physical s-p-coords
   * \param zeta linear coord point
   * \return derivative of itf at zeta
   */
  template <typename T>
  Eigen::Vector2<T> ditf(const Eigen::Ref<const Eigen::Vector2<T>> &zeta) const {
    return Eigen::Vector2<T>(0.5 / (std::sqrt(std::max((T)0., (T)1. - 1. / zeta(0))) * zeta(0) * zeta(0)),
                             p00 * ikappa * std::pow(std::max((T)0., (T)zeta(1)), ikappa - 1));
  }

  /*!
   * Transform a point physical s-p-coords -> linearized coords
   * \param x point in s-p-coords
   * \return x transformed to linearized coords
   */
  template <typename T>
  Eigen::Vector2<T> tf(const Eigen::Ref<const Eigen::Vector2<T>> &x) const {
    return Eigen::Vector2<T>(1. / (1. - x[0] * x[0]),
                             std::pow(x[1] / p00, kappa));
  }
  
  /*!
   * Derivative (diagonal of Jacobian, off-diagonal is 0) of transform physical s-p-coords -> linearized coords
   * \param x point in s-p-coords
   * \return derivative of tf at x
   */
  template <typename T>
  Eigen::Vector2<T> dtf(const Eigen::Ref<const Eigen::Vector2<T>> &x) const {
    return Eigen::Vector2<T>(2. * x[0] / ((1. - x[0] * x[0]) * (1. - x[0] * x[0])),
                             kappa / std::pow(p00, kappa) * std::pow(x[1], kappa - 1));
  }

  /*!
   * Return string representation (esp. for python)
   */
  std::string repr() const {
    std::stringstream ss;
    ss << "a = " << a
       << ", Omega = " << Omega
       << ", p00 = " << p00
       << ", kappa = " << kappa
       << ", cp = " << cp;
    return ss.str();
  }

  /*!
   * Convert into a python tuple for pickling
   * \param s physical_parameters instance to pickle
   * \return py::tuple containing relevant values of s
   */
  static py::tuple pickle(const physical_parameters &s) {
    return py::make_tuple(s.a, s.Omega, s.p00, s.kappa, s.cp);
  }

  /*!
   * Construct physical_parameters from python tuple for unpickling
   * \param t py::tuple containing values defining physical_parameters
   * \return physical_parameters defined by given tuple
   */
  static physical_parameters unpickle(py::tuple t) {
    if (t.size() != 5)
      throw std::runtime_error("invalid state in physical_parameters::unpickle");

    return physical_parameters(t[0].cast<double>(),
                               t[1].cast<double>(),
                               t[2].cast<double>(),
                               t[3].cast<double>(),
                               t[4].cast<double>());
  }
  
  /*!
   * Bind class to pybind11 module
   * \param m module to bind to
   */
  static void bind(py::module_ &m) {
    py::class_<physical_parameters>(m, "PhysicalParameters")
      .def(py::init<double, double, double, double, double>(),
           py::arg("a") = 6371000., py::arg("Omega") = 7.2921e-05,
           py::arg("p00") = 101325., py::arg("kappa") = 2. / 7., py::arg("cp") = 1003.5)
      .def(py::pickle(&physical_parameters::pickle, &physical_parameters::unpickle))
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
