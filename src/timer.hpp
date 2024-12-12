#ifndef TIMER_HPP
#define TIMER_HPP

#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>
#include "common.hpp"

class alignas(64) timer : public std::enable_shared_from_this<timer> {
public:
  std::vector<double> times; //!< list of accumulated times so far
  std::vector<std::string> names; //!< list of section names
  int current_section = -1; //!< currently active section or -1 if no section is active
  std::chrono::time_point<std::chrono::high_resolution_clock> t_last; //!< start of current measurement

  /*! Construct without any sections. */
  timer() {}

  /*! Construct with a given list of sections. */
  timer(std::vector<std::string> names) : times(names.size(), 0), names(names) { }

  /*! Reset all section times to 0. */
  void reset() {
    std::fill(times.begin(), times.end(), 0);
  }

  /*! Turn a string into a section index, if the section does not yet exist, create it.
   * \param s name of the section
   * \return segment index (to be used in start_section)
   */
  int get_index_from_name(std::string s) {
    int i = std::find(names.begin(), names.end(), s) - names.begin();
    if (i == (int)names.size())
      names.push_back(s), times.push_back(0);
    return i;
  }

  /*! Start timing section with index i, if a section is already active it will be ended at the same point in time.
   * \param i section index
   */
  void start_section(int i) {
    auto t = std::chrono::high_resolution_clock::now();
        
    double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
    if (current_section >= 0)
      times[current_section] += dt;
        
    t_last = t;
    current_section = i;
  }

  /*! End current section. */
  void end_section() {
    auto t = std::chrono::high_resolution_clock::now();
    double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
    if (current_section >= 0)
      times[current_section] += dt;
    current_section = -1;
  }

  /*! Format stored times with names into a human-readable string.
   * \param sep string to use to separate different sections (e.g. "\n" or ", ")
   * \return string containing measured times
   */
  std::string format_times(std::string sep = ", ") {
    auto t = std::chrono::high_resolution_clock::now();
    double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
    if (current_section >= 0)
      times[current_section] += dt;
    
    std::stringstream ss;
    ss << names[0] << ": " << times[0] * 1e-6 << "ms";
    for (int i = 1; i < (int)times.size(); i++)
      ss << sep << names[i] << ": " << times[i] * 1e-6 << "ms";
    return ss.str();
  }

  /*! Add times from another timer to this, matching sections of the same name and adding new section as necessary.
   * \param o other timer to add from
   * \return *this
   */
  timer &operator+=(const timer &to) {
    for (int i = 0; i < (int)to.names.size(); ++i) {
      int j = get_index_from_name(to.names[i]);
      times[j] += to.times[i];
    }
    return *this;
  }
  
  /*! Add two timers
    \sa timer::operator+=
   */
  friend timer operator+(timer lhs, const timer &rhs) {
    return lhs += rhs;
  }

  /*! Pack values into a python tuple for pickling support, note that current section is lost. */ 
   static py::tuple pickle(std::shared_ptr<timer> s) {
     return py::make_tuple(s->times, s->names);
   }

  /*! Unpack values from a python tuple for pickling support. */ 
   static std::shared_ptr<timer> unpickle(py::tuple t) {
     if (t.size() != 2)
       throw std::runtime_error("invalid state in timer::unpickle");
     
     auto res = std::make_shared<timer>(t[1].cast<std::vector<std::string>>());
     res->times = t[0].cast<std::vector<double>>();
     return res;
  }
 
  /*! Generate python bindings to given module. */
  static void bind(py::module_ &m) {
    py::class_<timer, std::shared_ptr<timer>>(m, "Timer")
      .def(py::init<std::vector<std::string>>(), py::arg("names"))
      .def_readonly("names", &timer::names)
      .def_readonly("times", &timer::times)
      .def("get_index_from_name", &timer::get_index_from_name)
      .def("start_section", &timer::start_section)
      .def("end_section", &timer::end_section)
      .def("format_times", &timer::format_times)
      .def(py::pickle(&timer::pickle, &timer::unpickle))
      .def(py::self + py::self)
      .def(py::self += py::self)
      .def("__repr__", [](std::shared_ptr<timer> t){ return t->format_times("\n"); });
  }
};

#endif
