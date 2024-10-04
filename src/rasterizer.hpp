#ifndef RASTERIZER_HPP
#define RASTERIZER_HPP

#include "common.hpp"
#include "pool.hpp"
#include "gfx/timsort.hpp"

class rasterizer
{
public:
  // calculate y-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at x-coord x
  static double y_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &x);
  // calculate x-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at y-coord y
  static double x_at(const double &x0, const double &y0,
                            const double &x1, const double &y1,
                            const double &y);
  // calculate the area of the intersection of the rectangle [xlo, xhi] x [ylo, yhi] with the
  // halfspace below the line connecting (xlo, syl) and (xhi, syr)
  static double rect_intersect_area(const double &xlo, const double &xhi,
                                           const double &ylo, const double &yhi,
                                           const double &syl, const double &syr);

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
            const int &idx, const int &sidx);    
    std::string repr(const int precision = -1) const;

    double slope() const;
    
    double y_at(const double &x) const;
    double x_at(const double &y) const;
    
    // check if there is a unique (non-endpoint-) intersect with o, if intersecting set p to intersection point and return true
    // p(0) =: x will be the largest value such that y_at(x) <=> o.y_at(x) is the same as at small x
    // p(1) =: y will be the average of y_at of the two segments
    bool intersect(segment &o, Eigen::Vector2d &p);

    static void bind(py::module_ &m);
  };

  struct segment_cmp {
    double *cur_x; // x-coordinate to use for comparison

    segment_cmp() = delete;
    segment_cmp(const segment_cmp &c);
    segment_cmp(double *cur_x);
    
    bool operator()(const segment &a, const segment &b) const;
    bool operator()(const segment *a, const segment *b) const;
  };

  struct event {
    void *p;

    event();
    event(segment &s, bool start);
    
    // get the position of the event
    const Eigen::Vector2d& pos() const;
    // get whether the event represents a segment start
    bool start() const;
    // get the segment
    segment* seg() const;

    friend bool operator<(const event &lhs, const event &rhs);
    std::string repr(const int &precision = -1) const;

    static void bind(py::module_ &m);
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

  void setup_timer();
#endif
  
  rasterizer() = default;

  rasterizer(const rasterizer &o) = delete;
  
  rasterizer(rasterizer &&o);
  
  // seg = (geometric) segments; seg[i] = {start point, end point, index of cell below segment}; segments should include top and bottom but not left and right border (direction does not matter)
  // bounds = (left, right, bottom, top) boundaries left < right bounds coordinate 0; bottom < top bounds coordinate 1
  rasterizer(const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &seg,
             const Eigen::Ref<const Eigen::Array4d> &bounds);

  ~rasterizer();

  // val = values of cells
  // res = resolution the result shall have
  Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
  rasterize(const Eigen::Ref<const Eigen::VectorXd> &val,
            const Eigen::Ref<const Eigen::Array2i> &res);

  static void bind(py::module_ &m);
protected:

  // find all intersects, split segments to remove intersects, update events accordingly
  // currently assumes there are no 3 segments intersecting in 1 point TODO: handle that case too
  int fix_intersects();
  
  // for given x, find pixel column containing it
  int x_to_pixel(const double &x) const;
  // for given y, gind pixel row containing it
  int y_to_pixel(const double &y) const;
  // for given i, give start of pixel column i
  double pixel_to_x(const int &i) const;
  // for given j, give start of pixel row j
  double pixel_to_y(const int &j) const;

  void add_area_between(const segment *s0, segment *s1, const double &x1);
};

#endif
