#ifndef RASTERIZER_HPP
#define RASTERIZER_HPP

#include "common.hpp"

class rasterizer
{
public:

  struct segment {
    Eigen::Vector2d s, e; // start and end point
    int idx, sidx;        // index of cell below (if it exists), index of segment
    double x0;            // (current) start of trapezoid below segment

    std::string repr(const int precision = -1) const {
      std::stringstream ss;
      if (precision > 0)
        ss << std::setprecision(precision);
      ss << "segment("
         << "s = (" << s(0) << ", " << s(1) << ")"
         << ", e = (" << e(0) << ", " << e(1) << ")"
         << ", idx = " << idx
         << ", sidx = " << sidx
         << ", x0 = " << x0
         << ")";
      return ss.str();
    }
  };

  typedef std::list<segment>::iterator sit_t;

  struct event {
    Eigen::Vector2d pos; // position of the event
    int n_in, n_out;     // number of incoming and outgoing segments
    sit_t its[3];        // segments associated to event (n_in incoming ones followed by n_out outgoing ones)
    friend bool operator<(const event &lhs, const event &rhs) {
      return std::tie(lhs.pos(0), lhs.pos(1)) < std::tie(rhs.pos(0), rhs.pos(1));
    }
    std::string repr(const int &precision = -1) const {
      std::stringstream ss;
      if (precision > 0)
        ss << std::setprecision(precision);
      ss << "event("
         << "coord = (" << pos(0) << ", " << pos(1) << ")"
         << ", in = [";
      for (int i = 0; i < n_in; ++i)
        ss << its[i]->sidx << ", ";
      ss << "]"
         << ", out = [";
      for (int i = 0; i < n_out; ++i)
        ss << its[i + n_in]->sidx << ", ";
      ss << "]"
         << ")";
      return ss.str();
    }
  };

  // segments that are in the current sweep line and ones that are not correspondingly
  std::list<segment> line, pool;
  std::vector<bool> used;

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
    used(seg.size()),
    bounds(bounds),
    lb(bounds(0), bounds(2)) {

    // costruct segments from input parameters and store them in pool
    std::vector<sit_t> sits;
    for (int i = 0; i < (int)seg.size(); ++i) {
      const auto &s = seg[i];
      if (std::get<0>(s)(0) < std::get<1>(s)(0)) {
        sits.push_back(pool.insert(pool.end(), {std::get<0>(s), std::get<1>(s), std::get<2>(s), i, std::get<0>(s)(0)}));
        // py::print("not reversed segment: ", *sits.back());
      } else if (std::get<0>(s)(0) > std::get<1>(s)(0)) {
        sits.push_back(pool.insert(pool.end(), {std::get<1>(s), std::get<0>(s), std::get<2>(s), i, std::get<1>(s)(0)}));
        // py::print("reversed segment :", *sits.back());
      } else {

        if (std::get<0>(s)(1) < std::get<1>(s)(1)) {
          sits.push_back(pool.insert(pool.end(), {std::get<0>(s), std::get<1>(s), std::get<2>(s), i, std::get<0>(s)(0)}));
        } else if (std::get<0>(s)(1) > std::get<1>(s)(1)) {
          sits.push_back(pool.insert(pool.end(), {std::get<1>(s), std::get<0>(s), std::get<2>(s), i, std::get<1>(s)(0)}));
        } else {
          throw std::runtime_error("rasterizer got 0-length segment, which is currently not supported");
        }

        // py::print(FORMAT("WARNING: rasterizer got exactly vertical segment {}", sits.back()->repr(16)));
      }
    }

#ifdef DEBUG_CHECKS
    for (const auto &it : sits)
      if (it->s(0) > it->e(0))
        throw std::runtime_error(FORMAT("segment was not oriented correctly : {}", it->repr(16)));
#endif

    // sort starting segments
    for (auto &i : start)
      pline.push_back(sits[i]);
    std::sort(pline.begin(), pline.end(), [](const sit_t &it0, const sit_t &it1){ return it0->s(1) < it1->s(1); });

    // collect and sort events
    for (auto &c : con) {
      if (c.size() < 2)
        throw std::runtime_error("rasterizer got invalid connection (size < 2)");
      if (c.size() > 3)
        throw std::runtime_error("rasterizer got invalid connection (size > 3)");

#ifdef DEBUG_CHECKS
      for (const int &i : c)
        if (i < 0 || i >= (int)sits.size())
          throw std::runtime_error(FORMAT("connection contains invalid index {}", i));
#endif

      // find common point
      const auto &s0 = sits[c[0]];
      const auto &s1 = sits[c[1]];
      Eigen::Vector2d p;
      if (s0->s == s1->s || s0->s == s1->e)
        p = s0->s;
      else if (s0->e == s1->s || s0->e == s1->e)
        p = s0->e;
      else
        throw std::runtime_error(FORMAT("rasterizer got invalid connection (no common point in first two segments {} and {})", s0->repr(), s1->repr()));

      // separate into incoming and outgoing segments
      std::vector<sit_t> seg_in, seg_out;
      for (const int &i : c) {
        if (sits[i]->s == p)
          seg_out.push_back(sits[i]);
        else if (sits[i]->e == p)
          seg_in.push_back(sits[i]);
        else
          throw std::runtime_error(FORMAT("rasterizer got invalid connection of size {}"
                                          "(participating segment {} does not contain common point ({}, {}))",
                                          con.size(), sits[i]->repr(), p(0), p(1)));
      }

#ifdef DEBUG_CHECKS
      if (seg_in.size() + seg_out.size() != c.size())
        throw std::runtime_error(FORMAT("incoming and outgoing segments don't add up: {} + {} != {}", seg_in.size(), seg_out.size(), c.size()));
#endif

      // sort segments by y-order (i.e. slope as they have a common end point)
      std::sort(seg_in.begin(), seg_in.end(),
                [](const sit_t &l, const sit_t &r){
                  return ((l->e(1) - l->s(1)) / (l->e(0) - l->s(0)) >
                          (r->e(1) - r->s(1)) / (r->e(0) - r->s(0))); });
      std::sort(seg_out.begin(), seg_out.end(),
                [](const sit_t &l, const sit_t &r){
                  return ((l->e(1) - l->s(1)) / (l->e(0) - l->s(0)) <
                          (r->e(1) - r->s(1)) / (r->e(0) - r->s(0))); });

      // #ifdef DEBUG_CHECKS
      //       for (int i = 0; i < (int)seg_in.size() - 1; ++i)
      //         if (!is_below(seg_in[i], seg_in[i+1])) {
      //           std::string msg = "segment ordering in event incorrect (seg_in): ";
      //           for (const auto &it : seg_in)
      //             msg += it->repr(16) + " ";
      //           throw std::runtime_error(msg);
      //         }
      //       for (int i = 0; i < (int)seg_out.size() - 1; ++i)
      //         if (!is_below(seg_out[i], seg_out[i+1])) {
      //           std::string msg = "segment ordering in event incorrect (seg_out): ";
      //           for (const auto &it : seg_out)
      //             msg += it->repr(16) + " ";
      //           throw std::runtime_error(msg);
      //         }
      // #endif

      // build event
      event e;
      e.pos = p;
      e.n_in = seg_in.size();
      e.n_out = seg_out.size();
      for (int i = 0; i < (int)seg_in.size(); ++i)
        e.its[i] = seg_in[i];
      for (int i = 0; i < (int)seg_out.size(); ++i)
        e.its[i + seg_in.size()] = seg_out[i];

      // add event to list
      events.push_back(e);
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

    // move all segments back into pool (will not do anything first time)
    std::fill(used.begin(), used.end(), false);
    pool.splice(pool.end(), line);
    // move start segments into line
    for (auto &it : pline) {
      it->x0 = lb(0);
      used[it->sidx] = true;
      line.splice(line.end(), pool, it);
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

protected:

  void process_event(const event &e) {
#ifdef DEBUG_CHECKS
    try {
      assert_line_order();
    } catch (std::runtime_error &ex) {
      throw std::runtime_error(FORMAT("before processing {}: {}", e.repr(15), ex.what()));
    }
#endif

    // add areas below ending segment
    for (int i = 0; i < e.n_in; ++i) {
      const auto &it = e.its[i];
      add_area_below(it, e.pos(0));
    }

    // find ins_pos (which is either the one following the highest ending segment in line or has to be found by linear search)
    sit_t ins_pos;
    if (e.n_in > 0) {
      ins_pos = std::next(e.its[e.n_in - 1]);
    } else {
      ins_pos = line.begin();
      while (ins_pos != line.end() && y_at(ins_pos, e.pos(0)) <= e.pos(1))
        ++ins_pos;
    }

    // add the area below ins_pos (i.e. area above topmost incoming segment or area which will be cut by event)
    if (ins_pos != line.end())
      add_area_below(ins_pos, e.pos(0));

    // remove ending segments from line
    for (int i = 0; i < e.n_in; ++i) {
      const auto &it = e.its[i];
#ifdef DEBUG_CHECKS
      if (!used[it->sidx])
        throw std::runtime_error(FORMAT("try to remove segment {} from line that is not contained", it->sidx));
#endif
      used[it->sidx] = false;
      pool.splice(pool.end(), line, it);
    }

    // add new segments to line (in order: starting at lowest, all before ins_pos)
    for (int i = 0; i < e.n_out; ++i) {
      const auto &it = e.its[i + e.n_in];
      it->x0 = e.pos(0);
#ifdef DEBUG_CHECKS
      if (used[it->sidx])
        throw std::runtime_error(FORMAT("try to add segment {} to line that is already contained", it->sidx));
#endif
      used[it->sidx] = true;
      line.splice(ins_pos, pool, it);
    }

#ifdef DEBUG_CHECKS
    try {
      assert_line_order();
    } catch (std::runtime_error &ex) {
      throw std::runtime_error(FORMAT("after processing (line length {}) {}: {}", line.size(), e.repr(15), ex.what()));
    }
#endif
  }

  // calculate y-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at x-coord x
  static inline double y_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &x) {
    if (x1 == x0 && x == x0) return 0.5 * (y1 + y0);
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }
  static inline double y_at(const sit_t &it, const double &x) {
    return y_at(it->s(0), it->s(1), it->e(0), it->e(1), x);
  }
  // calculate x-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at y-coord y
  static inline double x_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &y) {
    if (y1 == y0 && y == y0) return 0.5 * (x1 + x0);
    return x0 + (x1 - x0) * (y - y0) / (y1 - y0);
  }
  static inline double x_at(const sit_t &it, const double &y) {
    return x_at(it->s(0), it->s(1), it->e(0), it->e(1), y);
  }
  // check whether it0 is below it1 assuming they do not intersect except for possibly at endpoints and they have some common x
  static inline bool is_below(const sit_t &it0, const sit_t &it1) {
    if (it0->s(0) == it0->e(0)) {
      if (it1->s(0) == it1->e(0))
        return (it0->s(1) + it0->e(1)) < (it1->s(1) + it1->e(1));
      else
        return 0.5 * (it0->s(1) + it0->e(1)) < y_at(it1, it0->s(0));
    } else if (it1->s(0) == it1->e(0)) {
      return y_at(it0, it1->s(0)) < 0.5 * (it1->s(1) + it1->e(1));
    } else {
      double x = 0.5 * (std::max(it0->s(0), it1->s(0)) +  std::min(it0->e(0), it1->e(0)));
      return y_at(it0, x) < y_at(it1, x);
    }
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
    if (xlo >= xhi)
      return 0;
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

      if (xlo >= xhi)
        continue;

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
    if (line.size() <= 1)
      return;
    const double epsi_x = 1e-9, epsi_y = 1e-4;
    for (auto it1 = line.begin(), it0 = it1++; it1 != line.end(); it0 = it1++) {

      double sx = std::max(it0->s(0), it1->s(0));
      double ex = std::min(it0->e(0), it1->e(0));

      if (sx - ex > epsi_x)
        throw std::runtime_error(FORMAT("sweepline (length {}) contains segments without x-overlap: {} and {}",
                                        line.size(),
                                        it0->repr(16), it1->repr(16)));

      double x = 0.5 * (sx + ex);

      if (y_at(it0, x) - y_at(it1, x) > epsi_y)
        throw std::runtime_error(FORMAT("sweepline (length {}) is out of order between segments {} and {}, ys {} {} @ x {}",
                                        line.size(),
                                        it0->repr(16), it1->repr(16),
                                        y_at(it0, x), y_at(it1, x), x));
    }
  }
};

#define BIND_RASTERIZER_SEGMENT(m)                                      \
  py::class_<rasterizer::segment>(m, "RasterizerSegment")               \
  .def_readonly("s", &rasterizer::segment::s)                           \
  .def_readonly("e", &rasterizer::segment::e)                           \
  .def_readonly("idx", &rasterizer::segment::idx)                       \
  .def_readonly("sidx", &rasterizer::segment::sidx)                     \
  .def_readonly("x0", &rasterizer::segment::x0)                         \
  .def("__repr__", &rasterizer::segment::repr, py::arg("precision")=-1);


#define BIND_RASTERIZER_EVENT(m)                        \
  py::class_<rasterizer::event>(m, "RasterizerEvent")   \
  .def_readonly("pos", &rasterizer::event::pos)         \
  .def_readonly("n_in", &rasterizer::event::n_in)       \
  .def_readonly("n_out", &rasterizer::event::n_out)     \
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
  .def_readonly("pool", &rasterizer::pool)                              \
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
