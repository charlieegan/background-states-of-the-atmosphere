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

    double slope() const {
      return (e(1) - s(1)) / (e(0) - s(0));
    }
  };

  typedef std::list<segment>::iterator sit_t;

  struct event {
    //Eigen::Vector2d pos; // position of the event
    Eigen::Matrix<double, 2, 1> pos; // position of the event
    std::vector<sit_t> seg_in, seg_out; // incoming and outgoing segments
    friend bool operator<(const event &lhs, const event &rhs) {
      if (lhs.pos(0) != rhs.pos(0))
        return lhs.pos(0) < rhs.pos(0);

      for (auto &itl : lhs.seg_out)
        for (auto &itr : rhs.seg_in)
          if (itl == itr)
            return true;

      for (auto &itl : lhs.seg_in)
        for (auto &itr : rhs.seg_out)
          if (itl == itr)
            return false;

      return lhs.pos(1) < rhs.pos(1);
    }
    std::string repr(const int &precision = -1) const {
      std::stringstream ss;
      if (precision > 0)
        ss << std::setprecision(precision);
      ss << "event("
         << "coord = (" << pos(0) << ", " << pos(1) << ")"
         << ", in = [";
      for (auto &it : seg_in)
        ss << it->sidx << ", ";
      ss << "]"
         << ", out = [";
      for (auto &it : seg_out)
        ss << it->sidx << ", ";
      ss << "]"
         << ")";
      return ss.str();
    }
    void sort() {
      std::sort(seg_in.begin(), seg_in.end(),
                [](const sit_t &l, const sit_t &r){
                  return (l->slope() > r->slope()); });
      std::sort(seg_out.begin(), seg_out.end(),
                [](const sit_t &l, const sit_t &r){
                  return (l->slope() < r->slope()); });
    }
  };

  // segments that are in the current sweep line and ones that are not correspondingly
  std::list<segment> line, pool;
  std::vector<bool> used;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> out, fill;
  Eigen::Array4d bounds;
  Eigen::Array2d lb, step;
  double pixA;
  // pixel (i, j) represents area lb + (i, j)^T * step (pointwise multiplication)

  Eigen::VectorXd val; // value of function on cells

  std::vector<sit_t> pline; // starting line
  std::vector<event> events; // transition events

  rasterizer() = default;

  // seg = (geometric) segments; seg[i] = {start point, end point, index of cell below segment}; segments should include top and bottom but not left and right border (direction does not matter)
  // con = connectivity, con[i] has length >= 2 and contains indices of segments that are connected (order does not matter)
  // start = segment indices of segments at the start
  // bounds = (left, right, bottom, top) boundaries left < right bounds coordinate 0; bottom < top bounds coordinate 1
  rasterizer(const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &seg,
             const std::vector<std::vector<int>> &con,
             const std::vector<int> &start,
             const Eigen::Ref<const Eigen::Array4d> &bounds,
             const double &merge_epsi) :
    used(seg.size()),
    bounds(bounds),
    lb(bounds(0), bounds(2)) {

    // costruct segments from input parameters and store them in pool
    std::vector<sit_t> sits;
    for (int i = 0; i < (int)seg.size(); ++i) {
      const auto &s = seg[i];
      if (seg_is_oriented(std::get<0>(s), std::get<1>(s)))
        sits.push_back(pool.insert(pool.end(), {std::get<0>(s), std::get<1>(s), std::get<2>(s), i, std::get<0>(s)(0)}));
      else if (seg_is_oriented(std::get<1>(s), std::get<0>(s)))
        sits.push_back(pool.insert(pool.end(), {std::get<1>(s), std::get<0>(s), std::get<2>(s), i, std::get<1>(s)(0)}));
      else
        throw std::runtime_error("rasterizer got 0-length segment, which is currently not supported");
    }

#ifdef DEBUG_CHECKS
    for (const auto &it : sits)
      if (it->s(0) > it->e(0))
        throw std::runtime_error(FORMAT("segment was not oriented correctly : {}", it->repr(16)));
#endif

    // sort starting segments
    for (auto &i : start)
      pline.push_back(sits[i]);
    std::sort(pline.begin(), pline.end(), [](const sit_t &it0, const sit_t &it1){
      return std::make_pair(it0->s(1), it0->slope()) < std::make_pair(it1->s(1), it1->slope());
    });

    // collect events (based on connections)
    for (auto &c : con) {
      if (c.size() < 2)
        throw std::runtime_error("rasterizer got invalid connection (size < 2)");

#ifdef DEBUG_CHECKS
      for (const int &i : c)
        if (i < 0 || i >= (int)sits.size())
          throw std::runtime_error(FORMAT("connection contains invalid index {}", i));
#endif

      event e;

      // find common point
      const auto &s0 = sits[c[0]];
      const auto &s1 = sits[c[1]];
      if (s0->s == s1->s || s0->s == s1->e)
        e.pos = s0->s;
      else if (s0->e == s1->s || s0->e == s1->e)
        e.pos = s0->e;
      else
        throw std::runtime_error(FORMAT("rasterizer got invalid connection"
                                        "(no common point in first two segments {} and {})",
                                        s0->repr(), s1->repr()));

      // separate into incoming and outgoing segments
      for (const int &i : c) {
        if (sits[i]->s == e.pos)
          e.seg_out.push_back(sits[i]);
        else if (sits[i]->e == e.pos)
          e.seg_in.push_back(sits[i]);
        else
          throw std::runtime_error(FORMAT("rasterizer got invalid connection of size {}"
                                          "(participating segment {} does not contain common point ({}, {}))",
                                          con.size(), sits[i]->repr(), e.pos(0), e.pos(1)));
      }

#ifdef DEBUG_CHECKS
      if (e.seg_in.size() + e.seg_out.size() != c.size())
        throw std::runtime_error(FORMAT("incoming and outgoing segments don't add up: {} + {} != {}",
                                        e.seg_in.size(), e.seg_out.size(), c.size()));
#endif

      // add event to list
      events.push_back(e);
    }

    // sort events (by position / topology if neighboring)
    std::sort(events.begin(), events.end());

    // merge events which are within tolerance
    merge_close_events(merge_epsi);

    // sort event segments
    for (auto &e : events)
      e.sort();

    // find segment intersects (topology errors in input)
    auto ints = find_intersects();

    // fix topology errors by adding swap events at intersects
    if (!ints.empty()) {
      for (auto &e : ints)
        events.push_back(e);
      std::sort(events.begin(), events.end());
    }

#ifdef DEBUG_CHECKS
    if (!check_topology())
      throw std::runtime_error("invalid topology detected in rasterizer: " + top_error);
#endif
  }

  // val = values of cells
  // res = resolution the result shall have
  Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
  rasterize(const Eigen::Ref<const Eigen::VectorXd> &val,
            const Eigen::Ref<const Eigen::Array2i> &res) {

    this->val = val;

    // reset out
    out = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(res(0), res(1));
    fill = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(res(0), res(1));

    // calculate step spacings
    step = {(bounds(1) - bounds(0)) / res(0), (bounds(3) - bounds(2)) / res(1)};
    pixA = step(0) * step(1);

    // move all segments back into pool
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

    // add remaining areas at the end
    if (line.size())
      for (auto it1 = line.begin(), it0 = it1++; it1 != line.end(); it0 = it1++)
        add_area_between(it0, it1, bounds(1));

    return out;
  }

protected:

  std::vector<event> find_intersects() {
    std::vector<event> res;

    pool.splice(pool.end(), line);
    for (auto &it : pline) {
      it->x0 = lb(0);
      line.splice(line.end(), pool, it);
    }
    for (auto &e : events) {
      if (e.seg_in.empty() && e.seg_out.empty())
        continue;

      // insertion sort line and store inversions
      for (auto it1 = line.begin(), it0 = it1++; it1 != line.end(); ) {
        if (y_at(it0, e.pos(0)) > y_at(it1, e.pos(0))) {
          // found inversion
          // find intersection coordinate
          // (by binsearch to guarantee exact result: x is the first representable number at which the order breaks)
          double x = e.pos(0);
          for (double d = 1.0; x != (x - d); d *= 0.5)
            while (y_at(it0, x - d) > y_at(it1, x - d))
              x -= d;
          double y = y_at(std::abs(it0->slope()) < std::abs(it1->slope()) ? it0 : it1, x);
          // push event
          res.push_back({{x, y}, {it0, it1}, {it1, it0}});
          // swap segments in line
          line.splice(it0, line, it1);
          // if not at begin: go back a step to allow swaps over multiple positions
          if (it1 != line.begin())
            it0 = std::prev(it1);
          else
            it1 = std::next(it0);
        } else {
          it0 = it1++;
        }
      }
#ifdef DEBUG_CHECKS
      if (!std::is_sorted(line.begin(), line.end(),
                          [&e](const segment &lhs, const segment &rhs) {
                            return y_at(lhs, e.pos(0)) < y_at(rhs, e.pos(0));
                          }))
        throw std::runtime_error("wrong line order after insertion sort in find_intersects");
#endif

      sit_t ins_pos;
      if (e.seg_in.size()) {
        ins_pos = std::next(e.seg_in.back());
      } else {
        ins_pos = line.begin();
        while (ins_pos != line.end() && y_at(ins_pos, e.pos(0)) <= e.pos(1))
          ++ins_pos;
      }

      for (auto &it : e.seg_in) {
        pool.splice(pool.end(), line, it);
      }

      for (auto &it : e.seg_out) {
        it->x0 = e.pos(0);
        line.splice(ins_pos, pool, it);
      }

#ifdef DEBUG_CHECKS
      if (!std::is_sorted(line.begin(), line.end(),
                          [&e](const segment &lhs, const segment &rhs) {
                            return y_at(lhs, e.pos(0)) < y_at(rhs, e.pos(0));
                          }))
        throw std::runtime_error("wrong line order after event processing in find_intersects!");
#endif
    }

    return res;
  }

  std::string top_error;
  bool check_topology() {
    top_error = "";
    std::fill(used.begin(), used.end(), false);
    pool.splice(pool.end(), line);
    for (auto &it : pline) {
      it->x0 = lb(0);
      used[it->sidx] = true;
      line.splice(line.end(), pool, it);
    }
    for (int i = 0; i < (int)events.size(); ++i) {
      auto &e = events[i];
      
      if (e.seg_in.empty() && e.seg_out.empty())
        continue;

      sit_t ins_pos;
      if (e.seg_in.size()) {
        ins_pos = std::next(e.seg_in.back());
      } else {
        ins_pos = line.begin();
        while (ins_pos != line.end() && y_at(ins_pos, e.pos(0)) <= e.pos(1))
          ++ins_pos;
      }

      for (auto &it : e.seg_in) {
        if (!used[it->sidx]) {
          top_error = "double remove @ " + e.repr(16);
          return false;
        }
        used[it->sidx] = false;
        pool.splice(pool.end(), line, it);
      }

      for (auto &it : e.seg_out) {
        it->x0 = e.pos(0);
        if (used[it->sidx]) {
          top_error = "double insert @ " + e.repr(16);
          return false;
        }
        used[it->sidx] = true;
        line.splice(ins_pos, pool, it);
      }

      if ((i + 1 == (int)events.size() || events[i + 1].pos(0) != e.pos(0)) &&
          !std::is_sorted(line.begin(), line.end(),
                          [&e](const segment &lhs, const segment &rhs) {
                            return y_at(lhs, e.pos(0)) < y_at(rhs, e.pos(0));
                          })) {
        top_error = "order @ " + e.repr(16);
        return false;
      }
    }
    return true;
  }

  int merge_close_events(const double &merge_epsi) {
    int merge_cnt = 0;

    if (merge_epsi >= 0) {
      // merge very close events (and endpoints of corresponding segments)
      double epsi_x = merge_epsi * (bounds[1] - bounds[0]);
      double epsi_y = merge_epsi * (bounds[3] - bounds[2]);

      for (int l = 0, r = 0; r < (int)events.size(); ++l) {
        while (r < (int)events.size() && events[r].pos(0) - events[l].pos(0) <= epsi_x)
          ++r;

        // don't merge for empty events
        if (events[l].seg_in.empty() && events[l].seg_out.empty())
          continue;

        // find events to merge with
        for (int i = l + 1; i < r; ++i) {
          // check if event[i] is within merging distance
          if (std::abs(events[l].pos(1) - events[i].pos(1)) <= epsi_y) {
            // py::print("merge", events[l].repr(16), "and", events[i].repr(16));

            // merge events[i].seg_in; handle collapsed and flipped segments
            auto &p = events[l].pos;
            for (auto &it : events[i].seg_in) {

              it->e = p;

              // handle cases of segments collapsed to point or flipped due to merging
              if (seg_is_oriented(it->s,  it->e)) {
                // py::print("change endpoint of segment without flipping");
                events[l].seg_in.push_back(it);
              } else {
                int si = -1;
                for (int j = i - 1; j >= 0 && events[j].pos(0) >= it->s(0); --j) {
                  if (events[j].pos == it->s) {
                    si = j;
                    break;
                  }
                }
                if (si < 0)
                  throw std::runtime_error("did not find starting event of segment flipped in merging operation");

                std::erase(events[si].seg_out, it);

                if (it->s != it->e) {
                  // py::print("flipping segment");

                  std::swap(it->s, it->e);

                  events[l].seg_out.push_back(it);
                  events[si].seg_in.push_back(it);

                } else {
                  // py::print("segment collapsed to point is ignored");
                }
              }
            }
            events[i].seg_in.clear();

            // merge events[i].seg_out; segments here are always regular
            for (auto &it : events[i].seg_out) {
              it->s = p;
              events[l].seg_out.push_back(it);
            }
            events[i].seg_out.clear();

            ++merge_cnt;
          }
        }
      }
    }
    return merge_cnt;
  }

  void process_event(const event &e) {
    // ignore empty event
    if (e.seg_in.empty() && e.seg_out.empty())
      return;

    // add areas below ending segment
    for (auto &it : e.seg_in)
      add_area_below(it, e.pos(0));

    // find ins_pos (which is either the one following the highest ending segment in line or has to be found by linear search)
    sit_t ins_pos;
    if (e.seg_in.size()) {
      ins_pos = std::next(e.seg_in.back());
    } else {
      ins_pos = line.begin();
      while (ins_pos != line.end() && y_at(ins_pos, e.pos(0)) <= e.pos(1))
        ++ins_pos;
    }

    // add the area below ins_pos (i.e. area above topmost incoming segment or area which will be cut by event)
    if (ins_pos != line.end())
      add_area_below(ins_pos, e.pos(0));

    // remove ending segments from line
    for (auto &it : e.seg_in) {
#ifdef DEBUG_CHECKS
      if (!used[it->sidx])
        throw std::runtime_error(FORMAT("try to remove segment {} from line that is not contained", it->sidx));
#endif
      used[it->sidx] = false;
      pool.splice(pool.end(), line, it);
    }

    // add new segments to line (in order: starting at lowest, all before ins_pos)
    for (auto &it : e.seg_out) {
      it->x0 = e.pos(0);
#ifdef DEBUG_CHECKS
      if (used[it->sidx])
        throw std::runtime_error(FORMAT("try to add segment {} to line that is already contained", it->sidx));
#endif
      used[it->sidx] = true;
      line.splice(ins_pos, pool, it);
    }
  }

  // calculate y-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at x-coord x
  template <typename T = double>
  static inline double y_at(const T &x0, const T &y0,
                            const T &x1, const T &y1,
                            const T &x) {
    if (x1 == x0 && x == x0) return 0.5 * (y1 + y0);
    double res = y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    if (std::min(x0, x1) <= x && x <= std::max(x0, x1))
      res = std::clamp(res, std::min(y0, y1), std::max(y0, y1));
    return res;
  }
  template <typename T = double>
  static inline double y_at(const sit_t &it, const T &x) {
    return y_at<T>(*it, x);
  }
  template <typename T = double>
  static inline double y_at(const segment &s, const T &x) {
    return y_at<T>(s.s(0), s.s(1), s.e(0), s.e(1), x);
  }
  // calculate x-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at y-coord y
  static inline double x_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &y) {
    if (y1 == y0 && y == y0) return 0.5 * (x1 + x0);
    double res = x0 + (x1 - x0) * (y - y0) / (y1 - y0);
    if (std::min(y0, y1) <= y && y <= std::max(y0, y1))
      res = std::clamp(res, std::min(x0, x1), std::max(x0, x1));
    return res;
  }
  static inline double x_at(const sit_t &it, const double &y) {
    return x_at(*it, y);
  }
  static inline double x_at(const segment &s, const double &y) {
    return x_at(s.s(0), s.s(1), s.e(0), s.e(1), y);
  }

  static inline bool seg_is_oriented(const Eigen::Ref<const Eigen::Vector2d> &s, const Eigen::Ref<const Eigen::Vector2d> &e) {
    return (s(0) < e(0)) || ((s(0) == e(0)) && s(1) < e(1));
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

#ifdef DEBUG_CHECKS
    if (x0 > x1)
      throw std::runtime_error(FORMAT("add_area_between with negative x-interval {} to {}", x0, x1));
#endif

    // find pixel-index bounding box
    int imin = std::max(0,                   x_to_pixel(x0));
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
        fill(i, j) += std::clamp(A / pixA, 0.0, 1.0);
      }
    }

    // update current trapezoid starting x
    it1->x0 = x1;
  }

  void add_area_below(const sit_t &it, const double &x1) {
    // add area below it (if it exists)
    if (it == line.begin())
      return;
    add_area_between(std::prev(it), it, x1);
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


#define BIND_RASTERIZER_EVENT(m)                                        \
  py::class_<rasterizer::event>(m, "RasterizerEvent")                   \
  .def_readonly("pos", &rasterizer::event::pos)                         \
  .def_property_readonly("n_in", [](const rasterizer::event &e){ return e.seg_in.size(); }) \
  .def_property_readonly("n_out", [](const rasterizer::event &e){ return e.seg_out.size(); }) \
  .def("__repr__", &rasterizer::event::repr, py::arg("precision")=-1);


#define BIND_RASTERIZER(m)                                              \
  py::class_<rasterizer>(m, "Rasterizer")                               \
  .def(py::init<const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &, \
       const std::vector<std::vector<int>> &,                           \
       const std::vector<int> &,                                        \
       const Eigen::Ref<const Eigen::Array4d> &,                        \
       const double &>(),                                               \
       py::arg("seg"), py::arg("con"),                                  \
       py::arg("start"), py::arg("bounds"),                             \
       py::arg("merge_epsi"))                                           \
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
  .def_readonly("fill", &rasterizer::fill)                              \
  .def_readonly("bounds", &rasterizer::bounds)                          \
  .def_readonly("lb", &rasterizer::lb)                                  \
  .def_readonly("step", &rasterizer::step)                              \
  .def_readonly("pixA", &rasterizer::pixA);

#endif
