#ifndef RASTERIZER_HPP
#define RASTERIZER_HPP

#include "common.hpp"
#include "pool.hpp"
#include "gfx/timsort.hpp"

class rasterizer
{
public:
  // calculate y-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at x-coord x
  template <std::floating_point T>
  static inline T y_at(const T &x0, const T &y0,
                       const T &x1, const T &y1,
                       const T &x) {
// #ifdef DEBUG_CHECKS
//     if (x < x0 || x > x1) {
//       throw std::runtime_error(FORMAT("y_at evaluated outside segment ({} {}) -- ({} {}): {}\n{}",
//                                       x0, y0, x1, y1, x,
//                                       get_trace()));
//     }
// #endif
//    if (x0 == x1 && x == x0) return 0.5 * (y1 + y0);
    T t = (x - x0) / (x1 - x0);
    return (1 - t) * y0 + t * y1;
    // return std::lerp(y0, y1, t);
  }
  // calculate x-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at y-coord y
  template <std::floating_point T>
  static inline T x_at(const T &x0, const T &y0,
                       const T &x1, const T &y1,
                       const T &y) {
    if (y1 == y0 && y == y0) return 0.5 * (x1 + x0);
    T t = (y - y0) / (y1 - y0);
    return (1 - t) * x0 + t * x1;
    //return std::lerp(x0, x1, t);
  }

  struct segment_cmp;
  
  // alignment is needed for events
  struct alignas(64) segment {
    Eigen::Vector2d s, e;  // start and end point
    int idx, sidx;         // index of cell below (if it exists), index of segment
    double x0;             // (current) start of trapezoid below segment

    std::set<segment*>::iterator self; // set iterator (for line) of self (if it is in the line)

    segment() = default;
    
    segment(const Eigen::Ref<const Eigen::Vector2d> &s,
            const Eigen::Ref<const Eigen::Vector2d> &e,
            const int &idx, const int &sidx) :
      s(s), e(e),
      idx(idx), sidx(sidx),
      x0(s(0)) {
#ifdef DEBUG_CHECKS
      if (reinterpret_cast<long long>(this) & 0x3Fll)
        throw std::runtime_error(FORMAT("segment {} has incorrect alignment", sidx));
#endif
    }
    
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

    inline double slope() const {
      return (e(1) - s(1)) / (e(0) - s(0));
    }
    
    inline double y_at(const double &x) const {
      return rasterizer::y_at<double>(s(0), s(1), e(0), e(1), x);
    }
    inline double x_at(const double &y) const {
      return rasterizer::x_at<double>(s(0), s(1), e(0), e(1), y);
    }

    
    // check if there is a unique (non-endpoint-) intersect with o, if intersecting set p to intersection point and return true
    // p(0) =: x will be the largest value such that y_at(x) <=> o.y_at(x) is the same as at small x
    // p(1) =: y will be the average of y_at of the two segments
    bool inline intersect(segment &o, Eigen::Vector2d &p) {
      double xmin = std::max(s(0), o.s(0)), xmax = std::min(e(0), o.e(0));
      const double INF = std::numeric_limits<double>::infinity();

      // no x-intersection -> no intersect
      if (xmin >= xmax) {
        return false;
      }

      // y-coords at start and end
      double y0_t = y_at(xmin), y0_o = o.y_at(xmin);
      double y1_t = y_at(xmax), y1_o = o.y_at(xmax);

      auto ord0 = (y0_t <=> y0_o), ord1 = (y1_t <=> y1_o);
      
      // the segments do not intersect (case > 0) or touch at an endpoint (case == 0)
      if (ord0 * ord1 >= 0) {
        return false;
      }

      // TODO: use analytical solution to restrict search range for speedup
      // there is an intersect, find the approximate position approximately analytically
      // double a0 = y0_o - y0_t, a1 = y1_o - y1_t;
      // double x = std::lerp(xmin, xmax, a0 / (a0 - a1));

      // find the largest x s.t. (y_at(x) <=> o.y_at(x)) * ord0 > 0 by (hybrid) binary search starting at xmin,
      // i.e. the segments do not yet touch or have changed position
      double x = xmin, step = (xmax - xmin);
      for (; x + step != x; step *= 0.5)
        while (x + step <= xmax && (y_at(x + step) <=> o.y_at(x + step)) * ord0 > 0)
          x += step;

      // make sure we went all the way (make sure rounding at step addition did not screw us up)
      while ((y_at(std::nextafter(x, INF)) <=> o.y_at(std::nextafter(x, INF))) * ord0 > 0)
        x = std::nextafter(x, INF);
      
#ifdef DEBUG_CHECKS
      if (x < xmin || x > xmax)
        throw std::runtime_error("intersection is outside of x-bounds");

      if ((y_at(x) <=> o.y_at(x)) * ord0 < 0)
        throw std::runtime_error("calculated intersection is too far right (segments have swapped already)");
      
      if (std::nextafter(x, -INF) >= xmin && (y_at(std::nextafter(x, -INF)) <=> o.y_at(std::nextafter(x, -INF))) * ord0 <= 0)
        throw std::runtime_error("calculated intersection is too far right (at previous value segments already touch)");
#endif
      
      p(0) = x;
      p(1) = 0.5 * (y_at(x) + o.y_at(x));
      
      return true;
    }
  };

  struct segment_cmp {
    double *cur_x; // x-coordinate to use for comparison

    segment_cmp() = delete;
    segment_cmp(const segment_cmp &c) : cur_x(c.cur_x) {}
    segment_cmp(double *cur_x) : cur_x(cur_x) {}
    
    bool operator()(const segment &a, const segment &b) const {
      auto y0 = a.y_at(*cur_x), y1 = b.y_at(*cur_x);
      if (y0 != y1)
        return y0 < y1;
      if (a.slope() != b.slope()) {
        if (a.e == b.e)
          return a.slope() > b.slope();
        return a.slope() < b.slope();
      }
      return a.sidx < b.sidx;
    }
    bool operator()(const segment *a, const segment *b) const {
      return (*this)(*a, *b);
    }
  };

  struct event {
    void *p;

    event() : p(NULL) {}
    event(segment &s, bool start) : p(static_cast<void*>(start ? &s.s : &s.e)) {}
    
    // get the position of the event
    inline const Eigen::Vector2d& pos() const {
      return *static_cast<Eigen::Vector2d*>(p);
    }
    // get whether the event represents a segment start
    inline bool start() const {
      return !(reinterpret_cast<long long>(p) & 0x3Fll);
    }
    // get the segment
    inline segment* seg() const {
      return reinterpret_cast<segment*>(reinterpret_cast<long long>(p) & ~0x3Fll);
    }

    friend bool operator<(const event &lhs, const event &rhs) {
      // lower x first
      if (lhs.pos()(0) != rhs.pos()(0))
        return lhs.pos()(0) < rhs.pos()(0);

      // lower y first
      if (lhs.pos()(1) != rhs.pos()(1))
        return lhs.pos()(1) < rhs.pos()(1);
      
      // end events before start events
      if (lhs.start() != rhs.start())
        return lhs.start() < rhs.start();

      // lower slope first
      // if (lhs.seg()->slope() != rhs.seg()->slope())
      //   return lhs.seg()->slope() < rhs.seg()->slope();

      return lhs.seg()->sidx < rhs.seg()->sidx;
    }
    std::string repr(const int &precision = -1) const {
      std::stringstream ss;
      if (precision > 0)
        ss << std::setprecision(precision);
      ss << "event("
         << "coord = (" << pos()(0) << ", " << pos()(1) << ")"
         << (start() ? ", start segment " : ", end segment ")
         << seg()->sidx << ")";
      return ss.str();
    }
  };

  //std::list<segment> segs;                  // all segments (with begin left of end; vertical segments are left out)
  pool<segment, 1024> segs;                 // all segments (with begin left of end; vertical segments are left out)
  double *cur_x;                            // x-coord of sweep line
  segment_cmp cmp;                          // segment comparator
  std::set<segment*, segment_cmp> line;     // current sweep line
  std::vector<event> events;                // segment start / end events

  Eigen::VectorXd val;                      //  value of function on cells
  
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> out, fill;
  Eigen::Array4d bounds;
  Eigen::Array2d lb, step;
  double pixA;
  // pixel (i, j) represents area lb + (i, j)^T * step (pointwise multiplication)

#ifdef PROFILING
  std::shared_ptr<timer> time = NULL;
  int idx_time_prepare = -1;
  int idx_time_sort_events = -1;
  int idx_time_intersects_prepare = -1;
  int idx_time_intersects_neighbors = -1;
  int idx_time_intersects_intersect = -1;
  int idx_time_intersects_noint = -1;
  int idx_time_intersects_select = -1;
  int idx_time_intersects_erase_line = -1;
  int idx_time_intersects_erase_events = -1;
  int idx_time_intersects_make_segments = -1;
  int idx_time_intersects_upper_bound_iterate = -1;
  int idx_time_intersects_end_erase = -1;
  int idx_time_intersect_post = -1;
  int idx_time_rasterize_prepare = -1;
  int idx_time_rasterize_upper_bound = -1;
  int idx_time_rasterize_add_start = -1;
  int idx_time_rasterize_add_end = -1;
  int idx_time_rasterize_line_insert = -1;
  int idx_time_rasterize_line_erase = -1;

  void setup_timer() {
    if (!time)
      time = std::make_shared<timer>();
    idx_time_prepare = time->get_index_from_name("rasterizer.prepare");
    idx_time_sort_events = time->get_index_from_name("rasterizer.sort_events");
    idx_time_intersects_prepare = time->get_index_from_name("rasterizer.fix_intersects.prepare");
    idx_time_intersects_neighbors = time->get_index_from_name("rasterizer.fix_intersects.start.find_neighbors");
    idx_time_intersects_intersect = time->get_index_from_name("rasterizer.fix_intersects.start.calc_intersects");
    idx_time_intersects_noint = time->get_index_from_name("rasterizer.fix_intersects.start.no_intersect");
    idx_time_intersects_select = time->get_index_from_name("rasterizer.fix_intersects.start.select_intersect");
    idx_time_intersects_erase_line = time->get_index_from_name("rasterizer.fix_intersects.start.erase_line");
    idx_time_intersects_erase_events = time->get_index_from_name("rasterizer.fix_intersects.start.erase_events");
    idx_time_intersects_make_segments = time->get_index_from_name("rasterizer.fix_intersects.start.make_segments");
    idx_time_intersects_upper_bound_iterate = time->get_index_from_name("rasterizer.fix_intersects.start.upper_bound_iterate");
    idx_time_intersects_end_erase = time->get_index_from_name("rasterizer.fix_intersects.end.erase");
    idx_time_intersect_post = time->get_index_from_name("rasterizer.fix_intersects.post");
    idx_time_rasterize_prepare = time->get_index_from_name("rasterizer.rasterize.prepare");
    idx_time_rasterize_upper_bound = time->get_index_from_name("rasterizer.rasterize.upper_bound");
    idx_time_rasterize_add_start = time->get_index_from_name("rasterizer.rasterize.add_start");
    idx_time_rasterize_add_end = time->get_index_from_name("rasterizer.rasterize.add_end");
    idx_time_rasterize_line_insert = time->get_index_from_name("rasterizer.rasterize.insert");
    idx_time_rasterize_line_erase = time->get_index_from_name("rasterizer.rasterize.erase");
  }
#endif
  
  rasterizer() = default;

  rasterizer(const rasterizer &o) = delete;
  
  rasterizer(rasterizer &&o) :
    segs(std::move(o.segs)), cur_x(o.cur_x), cmp(cur_x), line(cmp),
    events(o.events),
    val(o.val), out(o.out), fill(o.fill),
    bounds(o.bounds), lb(o.lb), step(o.step)
#ifdef PROFILING
    , time(o.time), idx_time_prepare(o.idx_time_prepare)
#endif
  {
    //segs.splice(segs.begin(), o.segs);
    o.cur_x = NULL;
  }
  
  // seg = (geometric) segments; seg[i] = {start point, end point, index of cell below segment}; segments should include top and bottom but not left and right border (direction does not matter)
  // bounds = (left, right, bottom, top) boundaries left < right bounds coordinate 0; bottom < top bounds coordinate 1
  rasterizer(const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &seg,
             const Eigen::Ref<const Eigen::Array4d> &bounds) :
    cur_x(new double(bounds(0))),
    cmp(cur_x),
    line(cmp),
    bounds(bounds),
    lb(bounds(0), bounds(2)) {

#ifdef PROFILING
    setup_timer();
#endif

#ifdef PROFILING
    time->start_section(idx_time_prepare);
#endif
    
    // costruct segments and events from input parameters
    for (int i = 0; i < (int)seg.size(); ++i) {
      const auto &s = seg[i];
      if (std::get<0>(s)(0) < std::get<1>(s)(0)) {
        segs.push_back(segment(std::get<0>(s), std::get<1>(s), std::get<2>(s), (int)segs.size()));
        events.push_back(event(segs.back(), true));
        events.push_back(event(segs.back(), false));
      } else if (std::get<1>(s)(0) < std::get<0>(s)(0)) {
        segs.push_back(segment(std::get<1>(s), std::get<0>(s), std::get<2>(s), (int)segs.size()));
        events.push_back(event(segs.back(), true));
        events.push_back(event(segs.back(), false));
      } 
    }

#ifdef PROFILING
    time->start_section(idx_time_sort_events);
#endif
    
    // sort events
    //std::sort(events.begin(), events.end());
    gfx::timsort(events.begin(), events.end(), std::less<event>());
    
#ifdef PROFILING
    time->end_section();
#endif
    
    // fix intersects (by splitting segments)
    int fixcnt = fix_intersects();
    py::print(FORMAT("fixed {} intersects", fixcnt));

#ifdef EXPENSIVE_DEBUG_CHECKS
    int fixcnt2 = fix_intersects();
    py::print(FORMAT("second fix: {} intersects", fixcnt2));
#endif
  }

  ~rasterizer() {
    if (cur_x)
      delete cur_x;
  }

  // val = values of cells
  // res = resolution the result shall have
  Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
  rasterize(const Eigen::Ref<const Eigen::VectorXd> &val,
            const Eigen::Ref<const Eigen::Array2i> &res) {
    
#ifdef PROFILING
    setup_timer();
    time->start_section(idx_time_rasterize_prepare);
#endif
    
    this->val = val;
    
    *cur_x = bounds(0);
    // py::print(FORMAT("reset cur_x = {}", *cur_x));

    // reset out
    out = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(res(0), res(1));
    fill = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(res(0), res(1));

    // calculate step spacings
    step = {(bounds(1) - bounds(0)) / res(0), (bounds(3) - bounds(2)) / res(1)};
    pixA = step(0) * step(1);

#ifdef PROFILING
    time->end_section();
#endif

    
#ifdef DEBUG_CHECKS
    if (!line.empty())
      throw std::runtime_error("line not empty at start of rasterize");
#endif

    // process events in order
    for (auto &e : events) {
      *cur_x = e.pos()(0);
      
#ifdef EXPENSIVE_DEBUG_CHECKS
      for (auto &is : line)
        if (is->s(0) > *cur_x || is->e(0) < *cur_x)
          throw std::runtime_error(FORMAT("segment {} is in line that does not contain cur_x = {}", is->repr(), *cur_x));
      if (e.seg()->s(0) > *cur_x || e.seg()->e(0) < *cur_x)
        throw std::runtime_error(FORMAT("segment {} from event does not contain cur_x = {}", e.seg()->repr(), *cur_x));

      if (!std::is_sorted(line.begin(), line.end(), cmp)) {
        std::stringstream ss;
        for (auto &s : line) {
          ss << s->sidx << " ";
        }
        throw std::runtime_error(FORMAT("line {}is not sorted at x = {}, start of {}", ss.str(), *cur_x, e.repr()));
      }
#endif

#ifdef PROFILING
      time->start_section(idx_time_rasterize_upper_bound);
#endif

      // find segment above
      auto ub = e.start() ? line.upper_bound(e.seg()) : std::next(e.seg()->self);

#ifdef DEBUG_CHECKS
      if (!e.start() && ub == line.begin())
        throw std::runtime_error("ub == line.begin() at segment end");
#endif

#ifdef PROFILING
      time->start_section(e.start() ? idx_time_rasterize_add_start : idx_time_rasterize_add_end);
#endif

      // add areas that are topologically changed
      if (e.start()) {
        // when a segment starts, the area which it will split needs to be added (if there is such an area)
        if (ub != line.begin() && ub != line.end())
          add_area_between(*std::prev(ub), *ub, e.pos()(0));
      } else {
        // when a segment ends, the areas above and below it need to be added
        if (ub != line.end())
          add_area_between(e.seg(), *ub, e.pos()(0));
        if (e.seg()->self != line.begin())
          add_area_between(*std::prev(e.seg()->self), e.seg(), e.pos()(0));
      }

#ifdef PROFILING
      time->end_section();
#endif
      
      if (e.start()) {

#ifdef PROFILING
        time->start_section(idx_time_rasterize_line_insert);
#endif
        
        // if the event is a segment start, set starting x and insert segment
        e.seg()->x0 = e.pos()(0);
        e.seg()->self = line.insert(ub, e.seg());

#ifdef PROFILING
        time->end_section();
#endif
        
#ifdef DEBUG_CHECKS
        if (e.seg()->self != std::prev(ub))
          throw std::runtime_error("line insert at unexpected position");
#endif
      } else {
        // if the event is a segment end, std::prev(ub) will be the segment, erase it
#ifdef DEBUG_CHECKS
        if ((*std::prev(ub))->sidx != e.seg()->sidx)
          throw std::runtime_error("line erase at unexpected position");

        if (std::prev(ub) != line.begin() && e.seg()->x0 != e.pos()(0))
          throw std::runtime_error(FORMAT("segment not fully added when it is erased from line: x0 = {}, end = {}, missing {}",
                                          e.seg()->x0, e.pos()(0), e.pos()(0) - e.seg()->x0));
        
#endif

#ifdef PROFILING
        time->start_section(idx_time_rasterize_line_erase);
#endif

        line.erase(e.seg()->self);

#ifdef PROFILING
        time->end_section();
#endif
        
      }

#ifdef EXPENSIVE_DEBUG_CHECKS
      if (!std::is_sorted(line.begin(), line.end(), cmp)) {
        std::stringstream ss;
        for (auto &s : line) {
          ss << s->sidx << " ";
        }
        throw std::runtime_error(FORMAT("line {}is not sorted at x = {}, end of {}", ss.str(), *cur_x, e.repr()));
      }
#endif

    }

#ifdef DEBUG_CHECKS
    if (!line.empty())
      throw std::runtime_error("line not empty at the end");
#endif

    return out;
  }

protected:

  // find all intersects, split segments to remove intersects, update events accordingly
  // currently assumes there are no 3 segments intersecting in 1 point TODO: handle that case too
  int fix_intersects() {

#ifdef PROFILING
    time->start_section(idx_time_intersects_prepare);
#endif
    
    line.clear();
    *cur_x = bounds(0);

    std::set<event> aevents(events.begin(), events.end());

#ifdef PROFILING
    time->end_section();
#endif

#ifdef EXPENSIVE_DEBUG_CHECKS
    if (!std::is_sorted(aevents.begin(), aevents.end()))
      throw std::runtime_error("event set is not sorted, so comparison function must be messed up.");
#endif
    
    int icnt = 0;
    
    for (auto it = aevents.begin(); it != aevents.end(); ) {
      event e = *it;
      *cur_x = e.pos()(0);

      if (e.start()) {

#ifdef PROFILING
        time->start_section(idx_time_intersects_neighbors);
#endif

        // find neighboring segments in line
        segment *s = e.seg();
        auto lit = e.start() ? line.upper_bound(s) : std::next(s->self); // element that would come after s in line
        segment *s_prev = (lit == line.begin() ? NULL : *std::prev(lit));
        segment *s_next = (lit == line.end() ? NULL : *lit);

#ifdef PROFILING
        time->start_section(idx_time_intersects_intersect);
#endif
        
        // find possible intersects
        Eigen::Vector2d p_prev, p_next;
        bool int_prev = false, int_next = false;
        if (s_prev) int_prev = s->intersect(*s_prev, p_prev);
        if (s_next) int_next = s->intersect(*s_next, p_next);

#ifdef PROFILING
        time->end_section();
#endif
        
#ifdef DEBUG_CHECKS
        if (std::next(it) != aevents.upper_bound(e))
          throw std::runtime_error("upper bound is not next");
#endif
        
        if (!int_prev && !int_next) {
          // no intersects - just add segment to line and move on
#ifdef PROFILING
          time->start_section(idx_time_intersects_noint);
#endif
          s->self = line.insert(lit, s);
          ++it;
#ifdef PROFILING
          time->end_section();
#endif
        } else {
          // there is at least one intersect
          ++icnt;

#ifdef DEBUG_CHECKS
          if (int_next && int_prev && p_next(0) == p_prev(0))
            throw std::runtime_error("3 segments intersect in one point. This is not accounted for");
#endif

#ifdef PROFILING
          time->start_section(idx_time_intersects_select);
#endif
          
          // find which intersect happens first (if there are two)
          segment *s_int = NULL;
          Eigen::Vector2d p;
          if (int_next && (!int_prev || p_next(0) < p_prev(0))) {
            // next is the only or first intersect
            s_int = s_next;
            p = p_next;
          } else {
            // prev is the only or first intersect
            s_int = s_prev;
            p = p_prev;
          }

#ifdef DEBUG_CHECKS
          if (!s_int)
            throw std::runtime_error("failed to find intersecting segment");
#endif

#ifdef PROFILING
          time->start_section(idx_time_intersects_erase_line);
#endif          
          
          // delete intersecting segment from line
          line.erase(s_int->self);

#ifdef PROFILING
          time->start_section(idx_time_intersects_erase_events);
#endif          
          
          // delete old events
          aevents.erase(event(*s,     true )); // note: this is the event that is currently being processed, invalidating it
          aevents.erase(event(*s,     false)); // note: this event could be std::next(it), so just setting it to the
          aevents.erase(event(*s_int, true )); //       return value of the above erase will not do the trick
          aevents.erase(event(*s_int, false)); // 

#ifdef PROFILING
          time->start_section(idx_time_intersects_make_segments);
#endif          
          
          // create new split segments (if they are not vertical), add corresponding events,
          // add segment to line if the event should have already been processed:

          for (auto seg :
                 {segment(s->s, p, s->idx, 0),
                  segment(p, s->e, s->idx, 0),
                  segment(s_int->s, p, s_int->idx, 0),
                  segment(p, s_int->e, s_int->idx, 0)}) {

            if (seg.s(0) < seg.e(0)) {
              seg.sidx = (int)segs.size();
              segs.push_back(seg);
              aevents.insert(event(segs.back(), true));
              aevents.insert(event(segs.back(), false));
              if (event(segs.back(), true) < e)
                segs.back().self = line.insert(&segs.back()).first;
#ifdef DEBUG_CHECKS
              if (event(segs.back(), false) < e)
                throw std::runtime_error("end of split segment event is before current event");
#endif
            }
          }

#ifdef PROFILING
          time->start_section(idx_time_intersects_upper_bound_iterate);
#endif          
          // it was made invalid by deleting the old segment, so find the next position with upper_bound
          it = aevents.upper_bound(e);

#ifdef PROFILING
          time->end_section();
#endif          

#ifdef DEBUG_CHECKS
          if (it == aevents.end())
            throw std::runtime_error("did not find event after segment start");
#endif    
        }
      } else {
        // segment end event

#ifdef PROFILING
        time->start_section(idx_time_intersects_end_erase);
#endif
        
        line.erase(e.seg()->self);
        ++it;

#ifdef PROFILING
        time->end_section();
#endif
      }
    }

#ifdef DEBUG_CHECKS
    if (!line.empty()) {
      std::stringstream ss;
      for (auto &s : line) {
        ss << s->sidx << " ";
      }
      throw std::runtime_error(FORMAT("line not empty at the end : {}", ss.str()));
    }
#endif

#ifdef PROFILING
    time->start_section(idx_time_intersect_post);
#endif

    events = std::vector<event>(aevents.size());
    int i = 0;
    for (auto &e : aevents)
      events[i++] = e;
    //events = std::vector(aevents.begin(), aevents.end());

#ifdef PROFILING
    time->end_section();
#endif

    return icnt;
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

  void add_area_between(const segment *s0, segment *s1, const double &x1) {
    // add area between s0 and s1 (interpreted as lines) in range s0->x0 to x1 to all affected pixels
    // assume s0 is below s1 (in y-coords) and s1->x0 < x1
    // i.e. add val[s1->idx] to all pixels in the trapezoid (weighted towards 0 on the boundary)

    // py::print(FORMAT("add_area_between({}, {}, {})", s0->repr(), s1->repr(), x1));
    
    const double &x0 = s1->x0;

    if (x0 == x1)
      return;

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
      double y0l = s0->y_at(xlo), y0r = s0->y_at(xhi);
      double y1l = s1->y_at(xlo), y1r = s1->y_at(xhi);

      // find topology of intersections
      int jmin = std::max(0, y_to_pixel(std::min(y0l, y0r)));
      int jmax = std::min((int)out.cols() - 1, y_to_pixel(std::max(y1l, y1r)));

      for (int j = jmin; j <= jmax; ++j) {
        double yhi = pixel_to_y(j + 1), ylo = pixel_to_y(j);
        double A = (rect_intersect_area(xlo, xhi, ylo, yhi, y1l, y1r)
                    - rect_intersect_area(xlo, xhi, ylo, yhi, y0l, y0r));

#ifdef DEBUG_CHECKS
        if (s1->idx < 0 || s1->idx >= val.size())
          throw std::runtime_error(FORMAT("index in rasterizer out of bounds: {}", s1->idx));
#endif

        out(i, j) += std::clamp(A / pixA, 0.0, 1.0) * val[s1->idx];
        fill(i, j) += std::clamp(A / pixA, 0.0, 1.0);
      }
    }

    // update current trapezoid starting x
    s1->x0 = x1;
  }
};

#define BIND_RASTERIZER_SEGMENT(m)                                      \
  py::class_<rasterizer::segment>(m, "RasterizerSegment")               \
  .def_readonly("s", &rasterizer::segment::s)                           \
  .def_readonly("e", &rasterizer::segment::e)                           \
  .def_readonly("idx", &rasterizer::segment::idx)                       \
  .def_readonly("sidx", &rasterizer::segment::sidx)                     \
  .def_readonly("x0", &rasterizer::segment::x0)                         \
  .def("x_at", &rasterizer::segment::x_at, py::arg("y"))                \
  .def("y_at", &rasterizer::segment::y_at, py::arg("x"))                \
  .def("__repr__", &rasterizer::segment::repr, py::arg("precision")=-1);

#define BIND_RASTERIZER_EVENT(m)                                        \
  py::class_<rasterizer::event>(m, "RasterizerEvent")                   \
  .def_property_readonly("pos", &rasterizer::event::pos)                \
  .def_property_readonly("start", &rasterizer::event::start)            \
  .def_property_readonly("seg", [](const rasterizer::event &e){ return *e.seg(); }) \
  .def("__repr__", &rasterizer::event::repr, py::arg("precision")=-1);

#ifdef PROFILING
#define BIND_RASTERIZER_PROFILING()             \
  def_readonly("time", &rasterizer::time)
#else
#define BIND_RASTERIZER_PROFILING()
#endif

//   .def_readonly("segs", &rasterizer::segs)                            
#define BIND_RASTERIZER(m)                                              \
  py::class_<rasterizer>(m, "Rasterizer")                               \
  .def(py::init<const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &, \
                const Eigen::Ref<const Eigen::Array4d>&>(),             \
       py::arg("seg"), py::arg("bounds"))                               \
  .def("rasterize", &rasterizer::rasterize,                             \
       py::arg("val"), py::arg("res"))                                  \
  .def_property_readonly("line", [](const rasterizer &rast) {           \
    std::vector<rasterizer::segment> res;                               \
    for (auto &it : rast.line)                                          \
      res.push_back(*it);                                               \
    return res;                                                         \
  })                                                                    \
  .BIND_RASTERIZER_PROFILING()                                          \
  .def_readonly("events", &rasterizer::events)                          \
  .def_readonly("out", &rasterizer::out)                                \
  .def_readonly("fill", &rasterizer::fill)                              \
  .def_readonly("bounds", &rasterizer::bounds)                          \
  .def_readonly("lb", &rasterizer::lb)                                  \
  .def_readonly("step", &rasterizer::step)                              \
  .def_readonly("pixA", &rasterizer::pixA);

#endif
