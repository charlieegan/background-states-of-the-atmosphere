#ifndef RASTERIZER_HPP
#define RASTERIZER_HPP

#include "common.hpp"
#include "pool.hpp"
#include "gfx/timsort.hpp"


/*! A class to rasterize piecewise constant functions (e.g. on Laguerre diagrams) to a regular grid. */
class rasterizer
{
public:
  /*! Calculate y-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at x-coord x. */
  static double y_at(const double &x0, const double &y0,
                     const double &x1, const double &y1,
                     const double &x);
  /*! Calculate x-coordinate of segment (extended to line) (x0, y0) -- (x1, y1) at y-coord y. */
  static double x_at(const double &x0, const double &y0,
                     const double &x1, const double &y1,
                     const double &y);
  /*! Calculate the area of the intersection of the rectangle [xlo, xhi] x [ylo, yhi] with the
   *  halfspace below the line connecting (xlo, syl) and (xhi, syr)
   */
  static double rect_intersect_area(const double &xlo, const double &xhi,
                                    const double &ylo, const double &yhi,
                                    const double &syl, const double &syr);

  struct segment_cmp;
  
  /*! A segment between two points with auxillary information needed in the algorithm. */
  struct alignas(64) segment {
    Eigen::Vector2d s; //!< start point of segment
    Eigen::Vector2d e; //!< end point of segment (with e[1] > s[1])
    int idx;           //!< index of cell below segment
    int sidx;          //!< index of segment
    double x0;         //!< current start of trapezoid below segment
    std::set<segment*>::iterator self; //!< set iterator (for sweepline) of self (if it is in the sweepline)

    segment() = default;

    /*! Construct segment with given start, end, cell index and segment index. */
    segment(const Eigen::Ref<const Eigen::Vector2d> &s,
            const Eigen::Ref<const Eigen::Vector2d> &e,
            const int &idx, const int &sidx);

    /*! Return a string representation of the segment (mainly for python / debugging). */
    std::string repr(const int precision = -1) const;

    /*! Return the slope of the segment. */
    double slope() const;

    /*! Evaluate the y-coordinate of the segment at given x-coordinate. */
    double y_at(const double &x) const;

    /*! Evaluate the x-coordinate of the segment at given y-coordinate. */
    double x_at(const double &y) const;
    
    /*! Check if there is a unique (non-endpoint-) intersect with another segment o,
     *  if there is, set p to intersection point and return true.
     *  p(0) =: x will be the largest value such that y_at(x) <=> o.y_at(x) is the same as at small x
     *  p(1) =: y will be the average of y_at of the two segments
     */
    bool intersect(segment &o, Eigen::Vector2d &p);

    /*! Add python binding code to module m. */
    static void bind(py::module_ &m);
  };

  /*! A functor struct to compare two segments by their y-value at a point x. */
  struct segment_cmp {
    double *cur_x; //!< x-coordinate to use for comparison

    segment_cmp() = delete;
    segment_cmp(const segment_cmp &c);
    segment_cmp(double *cur_x);
    
    bool operator()(const segment &a, const segment &b) const;
    bool operator()(const segment *a, const segment *b) const;
  };

  /*! A struct to hold line-sweep events.
   *  Events are segment starts and segment ends, both are represented by a pointer
   *  to the corresponding point of the segment. Since segments are aligned by 64 bytes,
   *  if the event is a segment start, the pointer will have the 64 byte alignment,
   *  if the event is a segment end, the pointer will have an offset of sizeof(Eigen::Vector2d) == 16.
   */
  struct event {
    void *p; //*! pointer to start of end point of a segment

    event();
    
    /*! Construct start / end event from segment, based on value of start. */
    event(segment &s, bool start);
    
    /*! Get the position of the event. */
    const Eigen::Vector2d& pos() const;
    
    /*! Get whether the event represents a segment start. */
    bool start() const;
    
    /*! Get the segment the event is for. */
    segment* seg() const;

    /*! Compare two events by processing order (x-coord first, then y, then start/end, then index). */
    friend bool operator<(const event &lhs, const event &rhs);

    /*! Return a string representation of the event (mainly for python / debugging). */
    std::string repr(const int &precision = -1) const;

    /*! Add python bindings to module m. */
    static void bind(py::module_ &m);
  };

  pool<segment, 1024> segs;                 //!< all segments (with start.x < end.x; vertical segments are left out)
  double *cur_x;                            //!< current x-coord of sweep line
  segment_cmp cmp;                          //!< segment comparator
  std::set<segment*, segment_cmp> line;     //!< current sweep line
  std::vector<event> events;                //!< segment start / end events
  Eigen::VectorXd val;                      //!< value vector of function on cells
  
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> out;  //!< rasterized output
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> fill; //!< pixel fill amount in output
  Eigen::Array4d bounds;                    //!< domain boundaries [xmin, xmax, ymin, ymax]
  Eigen::Array2d lb;                        //!< left lower corner [xmin, ymin]
  Eigen::Array2d step;                      //!< pixel width and height
  double pixA;                              //!< pixel area
  // pixel (i, j) represents area lb + (i, j)^T * step (pointwise multiplication)

#ifdef PROFILING
  std::shared_ptr<timer> time = NULL;       //!< profiling timer
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
  
  /*! Construct rasterizer from segment list and bounds
   * \param seg (geometric) segments: seg[i] = {start point, end point, index of cell below segment}; segments should include top and bottom but do not have to include left and right border (direction of segments does not matter)
   * \param bounds rasterization boundaries (left, right, bottom, top) = (xmin, xmax, ymin, ymax)
   */
  rasterizer(const std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> &seg,
             const Eigen::Ref<const Eigen::Array4d> &bounds);

  ~rasterizer();

  /*! Rasterize a given function onto a regular grid with given resolution.
   *  This function can be called multiple times but the returned value (and set members) are only valid until the next execution (or destruction of the rasterizer).
   *  Make sure to take a copy if you intend to keep using the result beyond the lifetime of the rasterizer or when rasterizing multiple functions.
   * 
   * \param val value vector for the function; the indices are specified by the segments
   * \param resolution (xres, yres) that the result should have
   * 
   * \return matrix of shape (xres, yres) with the rasterized result. The result can also be retrieved by the member out and the member fill will contain the fill value of each pixel.
   */
  Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
  rasterize(const Eigen::Ref<const Eigen::VectorXd> &val,
            const Eigen::Ref<const Eigen::Array2i> &res);

  /*! Add python bindings to given module m. */
  static void bind(py::module_ &m);
protected:

  /*! Find all intersects, split segments to remove intersects, update events accordingly.
   *  Currently this functions assumes there are no 3 segments intersecting in a single point - TODO: handle that case too.
   * \return number of intersects fixed
   */
  int fix_intersects();
  
  /*! for given x, find pixel column containing it */
  int x_to_pixel(const double &x) const;
  
  /*! for given y, gind pixel row containing it */
  int y_to_pixel(const double &y) const;

  /*! for given i, give start of pixel column i */
  double pixel_to_x(const int &i) const;
  
  /*! for given j, give start of pixel row j */
  double pixel_to_y(const int &j) const;

  /*! Add up the are between two segments up to a given x-value to all affected pixels.
   * \param s0 lower segment
   * \param s1 upper segment (if s1 is not above s0 over the whole relevant range, behavior is undefined)
   * \param x1 x-value up to which to add up areas
   */
  void add_area_between(const segment *s0, segment *s1, const double &x1);
};

#endif
