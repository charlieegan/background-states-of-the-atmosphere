#ifndef TIMER_HPP
#define TIMER_HPP

#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>

class alignas(64) timer : public std::enable_shared_from_this<timer> {
public:
  std::vector<double> times;
  std::vector<std::string> names;
  int current_section = -1;
  std::chrono::time_point<std::chrono::high_resolution_clock> t_last;

  timer() {}
  timer(std::vector<std::string> names) : times(names.size(), 0), names(names) { }
  
  void reset() {
    std::fill(times.begin(), times.end(), 0);
  }

  int get_index_from_name(std::string s) {
    int i = std::find(names.begin(), names.end(), s) - names.begin();
    if (i == (int)names.size())
      names.push_back(s), times.push_back(0);
    return i;
  }
  
  void start_section(int i) {
    auto t = std::chrono::high_resolution_clock::now();
        
    double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
    if (current_section >= 0)
      times[current_section] += dt;
        
    t_last = t;
    current_section = i;
  }
    
  void end_section() {
    auto t = std::chrono::high_resolution_clock::now();
    double dt = std::chrono::duration_cast<std::chrono::nanoseconds>(t - t_last).count();
    if (current_section >= 0)
      times[current_section] += dt;
    current_section = -1;
  }
    
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


  timer &operator+=(const timer &to) {
    for (int i = 0; i < (int)to.names.size(); ++i) {
      int j = get_index_from_name(to.names[i]);
      times[j] += to.times[i];
    }
    return *this;
  }
  friend timer operator+(timer lhs, const timer &rhs) {
    return lhs += rhs;
  }
};

#define BIND_TIMER(m)                                                   \
  py::class_<timer, std::shared_ptr<timer>>(m, "Timer")                 \
  .def(py::init<std::vector<std::string>>(), py::arg("names"))          \
  .def_readonly("names", &timer::names)                                 \
  .def("get_index_from_name", &timer::get_index_from_name)              \
  .def("start_section", &timer::start_section)                          \
  .def("end_section", &timer::end_section)                              \
  .def("format_times", &timer::format_times)                            \
  .def(py::self + py::self)                                             \
  .def(py::self += py::self)                                            \
  .def("__repr__", [](std::shared_ptr<timer> t){ return t->format_times("\n"); });

#endif


