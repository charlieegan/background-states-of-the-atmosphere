#ifndef LAGUERRE_DIAGRAM_HPP
#define LAGUERRE_DIAGRAM_HPP

#include "common.hpp"

#include "physical_parameters.hpp"
#include "simulation_parameters.hpp"
#include "discretized_line_segment2.hpp"
#include "halfspace_intersection.hpp"
#include "rasterizer.hpp"

template <typename T> // dual type, used for duals and halfspace intersection (not everything)
class laguerre_diagram
{
public:
  typedef Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>> seeds_t;
  typedef Eigen::Vector3<T> Vec3;
  typedef Eigen::Vector4<T> Vec4;
  typedef Eigen::VectorX<T> VecX;


  struct diagram_edge {
    int pi, pj;
    int di, dj;

    discretized_line_segment<double> ls;

    // (approx) differential of area of pi (and negative area of pj) w.r.t. phi_pi - phi_pj
    double dif;

    std::string repr() const;

    static void bind(py::module_ &m);
  };

  struct virtual_diagram_edge {
    int pi, pj; // pi (with negative area) is connected to pj (with positive area)
    double dif; // differential
  };

  // input parameters
  int n;
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> ys;
  VecX duals;
  physical_parameters phys;
  simulation_parameters sim;

  // halfspace intersection
  halfspace_intersection<T> hs;
  std::vector<int> hints;
  
  // diagram:
  // vertices
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> verts;
  // edges
  std::vector<diagram_edge> edglist;
  // list of 'virtual' edges connecting negative areas with the closest positive ones
  std::vector<virtual_diagram_edge> vedglist;
  // primal faces (list of edges per face)
  std::vector<std::vector<int>> faces;
  // dual facets (list of edges per vertex)
  std::vector<std::vector<int>> dfacets;
  // face areas
  Eigen::VectorXd areas, areaerrs;

#ifdef PROFILING
  std::shared_ptr<timer> time = NULL;
  int idx_time_halfspace_intersection = -1,
    idx_time_extract_vertices = -1,
    idx_time_extract_edges = -1,
    idx_time_calc_areas = -1,
    idx_time_calc_integrals = -1,
    idx_time_prepare_rasterizer = -1,
    idx_time_jacobian = -1,
    idx_time_touching_dual = -1,
    idx_time_touching_dual_n = -1;

  void setup_timer();
#endif

  laguerre_diagram(const seeds_t &ys,
                   const Eigen::Ref<const VecX> &duals,
                   const physical_parameters &phys,
                   const simulation_parameters &sim);

  laguerre_diagram(const seeds_t &ys,
                   const Eigen::Ref<const VecX> &duals,
                   const physical_parameters &phys,
                   const simulation_parameters &sim,
                   const std::vector<int> &hints);

  
  void do_hs_intersect();

  void extract_diagram();

  VecX touching_dual(bool randomize);

  rasterizer get_rasterizer(std::function<Eigen::Vector2d(Eigen::Vector2d)> transform = identity());

  std::unordered_map<uint64_t, double> jac();

  std::pair<Eigen::VectorXd, std::pair<Eigen::VectorXi, Eigen::VectorXi>> jac_coo();

  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> get_poly(const int &pi_t);

  static void bind(py::module_ &m);

  static constexpr auto identity() {
    return std::function<Eigen::Vector2d(Eigen::Vector2d)>([](Eigen::Vector2d t){ return t; });
  }
};

#endif

