#ifndef RASTERIZER_HPP
#define RASTERIZER_HPP


#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <algorithm>

class rasterizer
{
public:

  struct segment {
    Eigen::Vector2d s, e; // start and end point
    int idx;              // index of cell below (if it exists)
    double x0;            // (current) start of trapezoid below segment

    std::string repr(const int precision = -1) const {
      std::stringstream ss;
      if (precision > 0)
        ss << std::setprecision(precision);
      ss << "segment("
         << "s = (" << s(0) << ", " << s(1) << ")"
         << ", e = (" << e(0) << ", " << e(1) << ")"
         << ", idx = " << idx
         << ", x0 = " << x0
         << ")";
      return ss.str();
    }
  };

  typedef std::list<segment>::iterator sit_t;

  enum event_type {
    TRANSITION, // segment 0 transitions to 1
    MERGE,      // segments 0, 1 merge into 2 (0 below 1)
    SPLIT,      // segment 0 splits into 1, 2 (1 below 2)
  };

  static std::string repr(const event_type &e) {
    switch (e) {
    case TRANSITION: return "TRANSITION";
    case MERGE: return "MERGE";
    case SPLIT: return "SPLIT";
    default: return "undefined";
    }
  }

  struct event {
    double pos;   // sweepline position (x-coord) of event
    event_type t; // type of event
    sit_t its[3]; // segments associated to event (for assignment see type definition)
    friend bool operator<(const event &lhs, const event &rhs) { return lhs.pos < rhs.pos; }
    std::string repr() const {
      std::stringstream ss;
      ss << "event("
         << "coord = (" << get_full_coord()(0) << ", " << get_full_coord()(1) << ")"
         << ", t = " << rasterizer::repr(t)
         << ")";
      return ss.str();
    }
    Eigen::Vector2d get_full_coord() const {
      return its[0]->e;
    }
  };

  // segments that were not yet visited, are completely visited and have intersection with the sweepline correspondingly
  std::list<segment> unused, used, line;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> out;
  Eigen::Array4d bounds;
  Eigen::Array2d lb, step;
  double pixA;
  // pixel (i, j) represents area lb + (i, j)^T * step (pointwise multiplication)

  Eigen::VectorXd val; // value of function on cells
  
  std::vector<sit_t> pline; // starting line
  std::vector<event> events; // transition events

  rasterizer() = default;
  
  // seg = (geometric) segments; seg[i] = {start point, end point, index of cell below segment}; segments should include top and bottom but not left and right border (direction does not matter)
  // con = connectivity, con[i] has length 2 or 3 and contains indices of segments that are connected (order does not matter)
  // start = segment indices of segments at the start
  // bounds = (left, right, bottom, top) boundaries left < right bounds coordinate 0; bottom < top bounds coordinate 1
  rasterizer(const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &seg,
             const std::vector<std::vector<int>> &con,
             const std::vector<int> &start,
             const Eigen::Ref<const Eigen::Array4d> &bounds) :
    bounds(bounds),
    lb(bounds(0), bounds(2)) {

    // costruct segments from input parameters and store them in unused
    std::vector<sit_t> sits;
    for (auto &s : seg) {
      if (std::get<0>(s)(0) < std::get<1>(s)(0)) {
        sits.push_back(unused.insert(unused.end(), {std::get<0>(s), std::get<1>(s), std::get<2>(s), std::get<0>(s)(0)}));
        // py::print("not reversed segment: ", *sits.back());
      } else if (std::get<0>(s)(0) > std::get<1>(s)(0)) {
        sits.push_back(unused.insert(unused.end(), {std::get<1>(s), std::get<0>(s), std::get<2>(s), std::get<1>(s)(0)}));
        // py::print("reversed segment :", *sits.back());
      } else {
        sits.push_back(unused.insert(unused.end(), {std::get<0>(s), std::get<1>(s), std::get<2>(s), std::get<0>(s)(0)}));
        throw std::runtime_error("rasterizer got exactly vertical segment");
      }
    }

    // move starting segments from unused to line (sorted by y-coord)
    for (auto &i : start)
      pline.push_back(sits[i]);
    std::sort(pline.begin(), pline.end(), [&](const sit_t &it0, const sit_t &it1){ return it0->s(1) < it1->s(1); });

    // collect and sort events
    for (auto &c : con) {
      if (c.size() == 2) {

        const auto &s0 = *sits[c[0]];
        const auto &s1 = *sits[c[1]];

        if (s0.e == s1.s) // end of segment 0 is start of segment 1
          events.push_back(event{s0.e(0), TRANSITION, {sits[c[0]], sits[c[1]], line.end()}});
        else if (s1.e == s0.s) // end of segment 1 is start of segment 0
          events.push_back(event{s1.e(0), TRANSITION, {sits[c[1]], sits[c[0]], line.end()}});
        else
          throw std::runtime_error("rasterizer got invalid transition");

      } else if (c.size() == 3) {

        const auto &s0 = sits[c[0]];
        const auto &s1 = sits[c[1]];
        const auto &s2 = sits[c[2]];

#define CHECK_COMMON(e0, e1, e2) ((e0 ? s0->e : s0->s) == (e1 ? s1->e : s1->s) && (e0 ? s0->e : s0->s) == (e2 ? s2->e : s2->s))

        if (CHECK_COMMON(1, 1, 0))
          if (is_below(s0, s1))
            events.push_back({s0->e(0), MERGE, {sits[c[0]], sits[c[1]], sits[c[2]]}});
          else
            events.push_back({s0->e(0), MERGE, {sits[c[1]], sits[c[0]], sits[c[2]]}});
        else if (CHECK_COMMON(1, 0, 1))
          if (is_below(s0, s2))
            events.push_back({s0->e(0), MERGE, {sits[c[0]], sits[c[2]], sits[c[1]]}});
          else
            events.push_back({s0->e(0), MERGE, {sits[c[2]], sits[c[0]], sits[c[1]]}});
        else if (CHECK_COMMON(0, 1, 1))
          if (is_below(s1, s2))
            events.push_back({s1->e(0), MERGE, {sits[c[1]], sits[c[2]], sits[c[0]]}});
          else
            events.push_back({s1->e(0), MERGE, {sits[c[2]], sits[c[1]], sits[c[0]]}});
        else if (CHECK_COMMON(1, 0, 0))
          if (is_below(s1, s2))
            events.push_back({s0->e(0), SPLIT, {sits[c[0]], sits[c[1]], sits[c[2]]}});
          else
            events.push_back({s0->e(0), SPLIT, {sits[c[0]], sits[c[2]], sits[c[1]]}});
        else if (CHECK_COMMON(0, 1, 0))
          if (is_below(s0, s2))
            events.push_back({s1->e(0), SPLIT, {sits[c[1]], sits[c[0]], sits[c[2]]}});
          else
            events.push_back({s1->e(0), SPLIT, {sits[c[1]], sits[c[2]], sits[c[0]]}});
        else if (CHECK_COMMON(0, 0, 1))
          if (is_below(s0, s1))
            events.push_back({s2->e(0), SPLIT, {sits[c[2]], sits[c[0]], sits[c[1]]}});
          else
            events.push_back({s2->e(0), SPLIT, {sits[c[2]], sits[c[1]], sits[c[0]]}});
        else
          throw std::runtime_error("rasterizer got invalid split/merge");

#undef CHECK_COMMON
#undef SEG_BELOW

      } else {
        throw std::runtime_error(FORMAT("rasterizer got connection with wrong number of arguments: {}", c.size()));
      }
    }
    std::sort(events.begin(), events.end());

  }

  // val = values of cells
  // res = resolution the result shall have
  Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
  rasterize(const Eigen::Ref<const Eigen::VectorXd> &val,
            const Eigen::Ref<const Eigen::Array2i> &res) {

    this->val = val;

    // reset out
    out = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(res(0), res(1));
    
    // calculate step spacings
    step = {(bounds(1) - bounds(0)) / res(0), (bounds(3) - bounds(2)) / res(1)};
    pixA = step(0) * step(1);

    // move all segments back into unused (will not do anything first time)
    unused.splice(unused.end(), line);
    unused.splice(unused.end(), used);
    // move start segments into line
    for (auto &it : pline) {
      it->x0 = lb(0);
      line.splice(line.end(), unused, it);
    }

    // process events in order
    for (auto &e : events)
      process_event(e);

    // add remaining areas for all left-over segments
    if (line.size())
      for (auto it1 = line.begin(), it0 = it1++; it1 != line.end(); it0 = it1++)
        add_area_between(it0, it1, bounds(1));

    return out;
  }
  
  void process_event(const event &e) {

    // py::print("process_event " + e.repr() + " line size: " + std::to_string(line.size()));

#ifdef DEBUG_CHECKS
    try {
      assert_line_order();
    } catch (std::runtime_error &ex) {
      throw std::runtime_error(FORMAT("before processing {}: {}", e.repr(), ex.what()));
    }
#endif
    
    // process line-sweep event e
    switch (e.t) {
    case TRANSITION: { // 0 -> 1
      add_area_below(e.its[0], e.pos);
      add_area_above(e.its[0], e.pos);
      e.its[1]->x0 = e.pos;
      line.splice(e.its[0], unused, e.its[1]);
      used.splice(used.end(), line, e.its[0]);
    } break;
    case MERGE: { // 0, 1 -> 2
      add_area_below(e.its[0], e.pos);
      add_area_between(e.its[0], e.its[1], e.pos);
      add_area_above(e.its[1], e.pos);
      e.its[2]->x0 = e.pos;
      line.splice(e.its[0], unused, e.its[2]);
      used.splice(used.end(), line, e.its[0]);
      used.splice(used.end(), line, e.its[1]);
    } break;
    case SPLIT: { // 0 -> 1, 2
      add_area_below(e.its[0], e.pos);
      add_area_above(e.its[0], e.pos);
      e.its[1]->x0 = e.pos;
      e.its[2]->x0 = e.pos;
      line.splice(e.its[0], unused, e.its[1]);
      line.splice(e.its[0], unused, e.its[2]);
      used.splice(used.end(), line, e.its[0]);
    } break;
    }

#ifdef DEBUG_CHECKS
    try {
      assert_line_order();
    } catch (std::runtime_error &ex) {
      throw std::runtime_error(FORMAT("after processing {}: {}", e.repr(), ex.what()));
    }
#endif
  }

  // calculate y-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at x-coord x
  static inline double y_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &x) {
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }
  static inline double y_at(const sit_t &it, const double &x) {
    return y_at(it->s(0), it->s(1), it->e(0), it->e(1), x);
  }
  // calculate x-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at y-coord y
  static inline double x_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &y) {
    return x0 + (x1 - x0) * (y - y0) / (y1 - y0);
  }
  static inline double x_at(const sit_t &it, const double &y) {
    return x_at(it->s(0), it->s(1), it->e(0), it->e(1), y);
  }
  // check whether it0 is below it1 assuming they do not intersect except for possibly at endpoints and they have some common x
  static inline bool is_below(const sit_t &it0, const sit_t &it1) {
    double x = 0.5 * (std::max(it0->s(0), it1->s(0)) +  std::min(it0->e(0), it1->e(0)));
    return y_at(it0, x) < y_at(it1, x);
  }
  // for given x, find pixel column containing it
  inline int x_to_pixel(const double &x) const {
    return std::floor((x - lb(0)) / step(0));
  }
  // for given y, gind pixel row containing it
  inline int y_to_pixel(const double &y) const {
    return std::floor((y - lb(1)) / step(1));
  }
  // for given i, give start of pixel column i
  inline double pixel_to_x(const int &i) const {
    return lb(0) + i * step(0);
  }
  // for given j, give start of pixel row j
  inline double pixel_to_y(const int &j) const {
    return lb(1) + j * step(1);
  }
  // calculate the area of the intersection of the rectangle [xlo, xhi] x [ylo, yhi] with the
  // halfspace below the line connecting (xlo, syl) and (xhi, syr)
  static inline double rect_intersect_area(const double &xlo, const double &xhi,
                                           const double &ylo, const double &yhi,
                                           const double &syl, const double &syr) {
    if (syl < ylo)
      if (syr < ylo) // both sides below -> no intersect
        return 0;
      else if (syr < yhi) // left side below, right inside -> lower right side triangle
        return 0.5 * (xhi - x_at(xlo, syl, xhi, syr, ylo)) * (syr - ylo);
      else // left side below, right above -> trapezoid with horizontal parallels on the right
        return 0.5 * ((xhi - x_at(xlo, syl, xhi, syr, ylo)) + (xhi - x_at(xlo, syl, xhi, syr, yhi))) * (yhi - ylo);
    else if (syl < yhi)
      if (syr < ylo) // left side inside, right below -> lower left side triangle
        return 0.5 * (x_at(xlo, syl, xhi, syr, ylo) - xlo) * (syl - ylo);
      else if (syr < yhi) // left side inside, right inside -> trepezoid with vertical parallels
        return 0.5 * (xhi - xlo) * ((syl - ylo) + (syr - ylo));
      else // left side inside, right above -> all - upper left triangle
        return (xhi - xlo) * (yhi - ylo) - 0.5 * (x_at(xlo, syl, xhi, syr, yhi) - xlo) * (yhi - syl);
    else
      if (syr < ylo) // left side above, right below -> trapezoid on left side with horizontal parallels
        return 0.5 * ((x_at(xlo, syl, xhi, syr, ylo) - xlo) + (x_at(xlo, syl, xhi, syr, yhi) - xlo)) * (yhi - ylo);
      else if (syr < yhi) // left side above, right inside -> all - upper right triangle
        return (xhi - xlo) * (yhi - ylo) - 0.5 * (xhi - x_at(xlo, syl, xhi, syr, yhi)) * (yhi - syr);
      else // left side above, right above -> completely contained
        return (xhi - xlo) * (yhi - ylo);
  }

  void add_area_between(const sit_t &it0, const sit_t &it1, const double &x1) {

    // add area between *it0 and *it1 (interpreted as lines) in range it1->x0 to x1 to all affected pixels
    // assume it0 is below it1 (in y-coords) and it1->x0 < x1
    // i.e. add val[it1->idx] to all pixels in the trapezoid (weighted towards 0 on the boundary)

    // py::print("add_area_between(" + it0->repr() + ", " + it1->repr() + ", " + std::to_string(x1) + ")");

    const double &x0 = it1->x0;

    // find pixel-index bounding box
    int imin = std::max(0,              x_to_pixel(x0));
    int imax = std::min((int)out.rows() - 1, x_to_pixel(x1));

    for (int i = imin; i <= imax; ++i) {
      // x left and right boundaries for pixel column
      double xlo = std::max(x0, pixel_to_x(i));
      double xhi = std::min(x1, pixel_to_x(i + 1));

      // y min-max range on left and right side
      double y0l = y_at(it0, xlo), y0r = y_at(it0, xhi);
      double y1l = y_at(it1, xlo), y1r = y_at(it1, xhi);

      // find topology of intersections
      int jmin = std::max(0, y_to_pixel(std::min(y0l, y0r)));
      int jmax = std::min((int)out.cols() - 1, y_to_pixel(std::max(y1l, y1r)));

      for (int j = jmin; j <= jmax; ++j) {
        double yhi = pixel_to_y(j + 1), ylo = pixel_to_y(j);
        double A = (rect_intersect_area(xlo, xhi, ylo, yhi, y1l, y1r)
                    - rect_intersect_area(xlo, xhi, ylo, yhi, y0l, y0r));

#ifdef DEBUG_CHECKS
        if (it1->idx < 0 || it1->idx >= val.size())
          throw std::runtime_error(FORMAT("index in rasterizer out of bounds: {}", it1->idx));
#endif
        
        out(i, j) += std::clamp(A / pixA, 0.0, 1.0) * val[it1->idx];
      }
    }

    // update current trapezoid starting x
    it1->x0 = x1;
  }

  void add_area_above(const sit_t &it, const double &x1) {
    // add area above it (if it exists)
    sit_t it1 = it; ++it1;
    if (it1 == line.end())
      return;
    add_area_between(it, it1, x1);
  }

  void add_area_below(const sit_t &it, const double &x1) {
    // add area below it (if it exists)
    if (it == line.begin())
      return;
    sit_t it0 = it; --it0;
    add_area_between(it0, it, x1);
  }

  void assert_line_order() {
    const double epsi = 1e-9;
    for (auto it1 = line.begin(), it0 = it1++; it1 != line.end(); it0 = it1++) {

      double sx = std::max(it0->s(0), it1->s(0));
      double ex = std::min(it0->e(0), it1->e(0));

      if (sx - ex > epsi)
        throw std::runtime_error(FORMAT("sweepline contains segments without x-overlap: {} and {}",
                                        it0->repr(16), it1->repr(16)));

      double x = 0.5 * (sx + ex);
      
      if (y_at(it0, x) - y_at(it1, x) > epsi)
        throw std::runtime_error(FORMAT("sweepline is out of order between segments {} and {}, ys {} {} @ x {}",
                                        it0->repr(16), it1->repr(16),
                                        y_at(it0, x), y_at(it1, x), x));
    }
  }
};

#define BIND_RASTERIZER_SEGMENT(m)                              \
  py::class_<rasterizer::segment>(m, "RasterizerSegment")       \
  .def_readonly("s", &rasterizer::segment::s)                   \
  .def_readonly("e", &rasterizer::segment::e)                   \
  .def_readonly("idx", &rasterizer::segment::idx)               \
  .def_readonly("x0", &rasterizer::segment::x0)                 \
  .def("__repr__", &rasterizer::segment::repr);

#define BIND_RASTERIZER_EVENT_TYPE(m)                           \
  py::enum_<rasterizer::event_type>(m, "RasterizerEventType")   \
  .value("TRANSITION", rasterizer::event_type::TRANSITION)      \
  .value("MERGE", rasterizer::event_type::MERGE)                \
  .value("SPLIT", rasterizer::event_type::SPLIT)                \
  .export_values();

#define BIND_RASTERIZER_EVENT(m)                                        \
  py::class_<rasterizer::event>(m, "RasterizerEvent")                   \
  .def_readonly("pos", &rasterizer::event::pos)                         \
  .def_readonly("t", &rasterizer::event::t)                             \
  .def("get_full_coord", &rasterizer::event::get_full_coord)            \
  .def("__repr__", &rasterizer::event::repr);


#define BIND_RASTERIZER(m)                                              \
  py::class_<rasterizer>(m, "Rasterizer")                               \
  .def(py::init<const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &, \
       const std::vector<std::vector<int>> &,                           \
       const std::vector<int> &,                                        \
       const Eigen::Ref<const Eigen::Array4d> &>(),                     \
       py::arg("seg"), py::arg("con"),                                  \
       py::arg("start"), py::arg("bounds"))                             \
  .def("rasterize", &rasterizer::rasterize,                             \
       py::arg("val"), py::arg("res"))                                  \
  .def_readonly("unused", &rasterizer::unused)                          \
  .def_readonly("used", &rasterizer::used)                              \
  .def_readonly("events", &rasterizer::events)                          \
  .def_property_readonly("pline", [](const rasterizer &rast) {          \
    std::vector<rasterizer::segment> res;                               \
    for (auto &it : rast.pline)                                         \
      res.push_back(*it);                                               \
    return res;                                                         \
  })                                                                    \
  .def_readonly("line", &rasterizer::line)                              \
  .def_readonly("out", &rasterizer::out)                                \
  .def_readonly("bounds", &rasterizer::bounds)                          \
  .def_readonly("lb", &rasterizer::lb)                                  \
  .def_readonly("step", &rasterizer::step)                              \
  .def_readonly("pixA", &rasterizer::pixA);

#endif
