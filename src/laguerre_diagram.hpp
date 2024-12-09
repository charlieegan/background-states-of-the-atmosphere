#ifndef LAGUERRE_DIAGRAM_HPP
#define LAGUERRE_DIAGRAM_HPP

#include "common.hpp"

#include "physical_parameters.hpp"
#include "simulation_parameters.hpp"
#include "discretized_line_segment.hpp"
#include "halfspace_intersection.hpp"
#include "rasterizer.hpp"

/*!
 * A class to calculate Laguerre diagrams.
 *
 * Given seed positions, dual values, physical parameters (defining the cost function) and
 * simulation parameters (defining among others the diagram boundaries and resolution), this class
 * allows calculating the Laguerre diagram including the geometry, cell areas and Jacobian of cell
 * areas w.r.t. the duals.
 */
template <typename T> // dual type, used for duals and halfspace intersection (not everything)
class laguerre_diagram : public std::enable_shared_from_this<laguerre_diagram<T>>
{
public:
  typedef Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>> seeds_t;
  typedef Eigen::Vector3<T> Vec3;
  typedef Eigen::Vector4<T> Vec4;
  typedef Eigen::VectorX<T> VecX;


  /*!
   * A struct representing a primal-dual edge with attached geometry.
   */
  struct diagram_edge {
    int pi; //!< primal staring index (first 6 are boundaries, next n are physical cells, last are top cells)
    int pj; //!< primal ending index (first 6 are boundaries, next n are physical cells, last are top cells)
    int di; //!< dual staring index 
    int dj; //!< dual ending index 

    discretized_line_segment<double> ls; //!< line geometry of (dual) edge in physical space

    double dif; //!< (approx) differential of area of primal cell with index pi (or pj with opposite sign) w.r.t. (phi_pi - phi_pj)

    /*!
     * Return string representation (esp. for python)
     */
    std::string repr() const;

    /*!
     * Bind class to pybind11 module
     * \param m module to bind to
     */
    static void bind(py::module_ &m);
  };

  /*!
   * A struct representing non-geometric edges for area modifications.
   */
  struct virtual_diagram_edge {
    int pi; //!< primal staring index (first 6 are boundaries, next n are physical cells, last are top cells)
    int pj; //!< primal ending index (first 6 are boundaries, next n are physical cells, last are top cells)
    double dif; //!< differential of area of primal cell with index pi (or pj with opposite sign) w.r.t. (phi_pi - phi_pj)
  };

  std::shared_ptr<laguerre_diagram<T>> parent; //!< parent diagram; when not NULL, negative areas will be kept consistent with parent

  int n; //!< number of seeds
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> ys; //!< seed positions
  VecX duals; //!< duals corresponding to seeds
  physical_parameters phys; //!< physical parameters defining cost / coordinate transform
  simulation_parameters sim; //!< simulation parameters defining boundaries, resolutions, etc.

  halfspace_intersection<T> hs; //!< halfspace intersection solver used to calculate diagram topology
  std::vector<int> hints; //!< hints for halfspace_intersection (primal indices close to the intersection for each halfspace)
  
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> verts; //!< geometric diagram vertices 
  std::vector<diagram_edge> edglist; //!< diagram edges (primal-dual with geometry)
  std::vector<virtual_diagram_edge> vedglist; //!< non-geometrical diagram edges (only primal with differential)
  std::vector<std::vector<int>> faces; //!< diagram faces (each face being represented by a list of edge indices)
  std::vector<std::vector<int>> dfacets; //!< dual facets, i.e. geometric vertices (each dual facet being represented by list of edge indices)
  Eigen::VectorXd areas; //!< diagram face areas
  Eigen::VectorXd areaerrs; //!< diagram face area error bounds (based on discretization)

#ifdef PROFILING
  std::shared_ptr<timer> time = NULL; //!< timer used for profiling
  int idx_time_halfspace_intersection = -1,
    idx_time_extract_vertices = -1,
    idx_time_extract_edges = -1,
    idx_time_calc_areas = -1,
    idx_time_calc_integrals = -1,
    idx_time_prepare_rasterizer = -1,
    idx_time_jacobian = -1,
    idx_time_touching_dual = -1,
    idx_time_touching_dual_n = -1;

  /*! Setup the profiling timer */
  void setup_timer();
#endif

  /*!
   * Construct Laguerre diagram based on seed positions, dual values, physical- and simulation parameters.
   * \param ys seed positions
   * \param duals dual values (including top cell, should be one longer than ys)
   * \param phys physical parameters, defining cost function / coordinate transform
   * \param sim simulation parameters, defining boundaries, resolution, etc.
   */
  laguerre_diagram(const seeds_t &ys,
                   const Eigen::Ref<const VecX> &duals,
                   const physical_parameters &phys,
                   const simulation_parameters &sim);

  /*!
   * Construct Laguerre diagram based on parent diagram with new duals
   * \param parent diagram to take values from
   * \param duals new dual values (should be of the same length as for parent)
   * \param keep_parent if true, store the parent; if false only take values from the parent but do not keep the link
   */
  laguerre_diagram(std::shared_ptr<laguerre_diagram<T>> parent,
                   const Eigen::Ref<const VecX> &duals,
                   const bool &keep_parent=true);

  /*! Detach from parent (set parent to NULL) */
  void detach();

  /*!
   * Calculate a dual vector where for each primal cell, if it was the only one changed, it is non-empty.
   *
   * For each primal index with an empty area, calculate the values of the corresponding dual at which the area becomes non-negative
   * and the value at which a neighboring area becomes 0. Return either the mean of these two values (if randomize == False) or
   * a random value between 25% and 75% of this interval.
   * \param randomize whether to return a random result
   */
  VecX touching_dual(bool randomize);

  /*!
   * Get a rasterizer corresponding to the diagram.
   * 
   * The rasterizer allows approximating functions defined piecewise constant over the diagram cells on a regular grid.
   * Optionally a transform can be given, which is applied to the points of the diagram before passing them to the rasterizer.
   * \param transform function to apply to points; this should act on the two dimensions independently and monotonously.
   * \return rasterizer instance corresponding to (transformed) diagram
   */
  rasterizer get_rasterizer(std::optional<std::function<Eigen::Vector2d(Eigen::Vector2d)>> transform = std::nullopt);

  /*!
   * Get the Jacobian of areas w.r.t. duals.
   *
   * This returns the Jacobian in coo-format as expected by scipy coo_array.
   * \return {v, {i, j}} such that J[i[k], j[k]] = v[k] represents the Jacobian
   */
  std::pair<Eigen::VectorXd, std::pair<Eigen::VectorXi, Eigen::VectorXi>> jac_coo();

  /*!
   * Get a polygon representing a diagram cell.
   *
   * \param pi_t index (corresponding to seeds / duals) of cell to get polygon for
   * \return n x 2 matrix with vertex coordinates
   */
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> get_poly(const int &pi_t);

  /*!
   * Bind class to pybind11 module
   * \param m module to bind to
   */
  static void bind(py::module_ &m);

private:
  
  /*!
   * Get the Jacobian of areas w.r.t. duals.
   * \return unordered_map from indices (combined into lower and upper word of a uint64) to values
   */
  std::unordered_map<uint64_t, double> jac();

  /*! Perform the halfspace intersection as part of diagram calculation */
  void do_hs_intersect();

  /*! Extract the diagram topology and geometry, areas and differentials from halfspace intersection */
  void extract_diagram();

};

#endif

