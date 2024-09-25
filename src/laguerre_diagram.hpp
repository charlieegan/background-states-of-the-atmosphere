#ifndef LAGUERRE_DIAGRAM_HPP
#define LAGUERRE_DIAGRAM_HPP

#include "common.hpp"

#include "physical_parameters.hpp"
#include "simulation_parameters.hpp"
#include "discretized_line_segment.hpp"
#include "halfspace_intersection.hpp"
#include "rasterizer.hpp"

template <typename T> // dual type, used for duals and halfspace intersection (not everything)
class laguerre_diagram
{
  static inline double sqr(const double &x) { return x * x; }

public:
  typedef Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>> seeds_t;
  typedef Eigen::Vector3<T> Vec3;
  typedef Eigen::Vector4<T> Vec4;
  typedef Eigen::VectorX<T> VecX;


  struct diagram_edge {
    int pi, pj;
    int di, dj;

    discretized_line_segment ls;

    // (approx) differential of area of pi (and negative area of pj) w.r.t. phi_pi - phi_pj
    double dif;

    std::string repr() const {
      std::stringstream ss;
      ss << "diagram_edge("
         << "pi = " << pi
         << ", pj = " << pj
         << ", di = " << di
         << ", dj = " << dj
         << ")";
      return ss.str();
    }
  };

  // input parameters
  int n;
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> ys;
  VecX duals;
  physical_parameters phys;
  simulation_parameters sim;

  // halfspace intersection
  halfspace_intersection<T> hs;

  // diagram:
  // vertices
  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> verts;
  // edges
  std::vector<diagram_edge> edglist;
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
    idx_time_jacobian = -1;

  void setup_timer() {
    if (!time)
      time = std::make_shared<timer>();
    idx_time_halfspace_intersection = time->get_index_from_name("laguerre_diagram.halfspace_intersection");
    idx_time_extract_vertices = time->get_index_from_name("laguerre_diagram.extract_vertices");
    idx_time_extract_edges = time->get_index_from_name("laguerre_diagram.extract_edges");
    idx_time_calc_areas = time->get_index_from_name("laguerre_diagram.calc_areas");
    idx_time_calc_integrals = time->get_index_from_name("laguerre_diagram.calc_integrals");
    idx_time_jacobian = time->get_index_from_name("laguerre_diagram.jacobian");
    idx_time_prepare_rasterizer = time->get_index_from_name("laguerre_diagram.prepare_rasterizer");
  }
#endif

  laguerre_diagram(const seeds_t &ys,
                   const Eigen::Ref<const VecX> &duals,
                   const physical_parameters &phys,
                   const simulation_parameters &sim) :
    n(ys.rows()), ys(ys), duals(duals), phys(phys), sim(sim) {

#ifdef PROFILING
    setup_timer();
#endif

    if ((int)duals.size() != n + 1)
      throw std::runtime_error("the size of the dual vector must be one larger than the first dimension of seed positions (the last entry of duals corresponds to the boundary cell)");

    do_hs_intersect();

    extract_diagram();
  }

  void do_hs_intersect() {
#ifdef PROFILING
    time->start_section(idx_time_halfspace_intersection);
#endif

    // perform the halfspace intersection

    Vec3 zeta_min, zeta_max;
    zeta_min << phys.tf<T>(sim.spmin.cast<T>()), (T)-1e100;
    zeta_max << phys.tf<T>(sim.spmax.cast<T>()), (T)1e100;
    hs = halfspace_intersection<T>(zeta_min, zeta_max);

    // add seed halfspaces
    for (int i = 0; i < n; i++) {
      Vec4 H(-0.5 * sqr(ys(i, 0) / phys.a),
             -phys.cp * ys(i, 1),
             1.0,
             phys.Omega * ys(i, 0) + duals(i));

      hs.add_halfspace(H);
    }

    // add boundary halfspaces
    for (int i = 0; i < sim.boundary_res; i++) {
      T s = sim.spmin(0) + (sim.spmax(0) - sim.spmin(0)) * i / (sim.boundary_res - 1);
      T z0 = 1. / (1 - s * s);

      Vec4 H(-0.5 * sqr(phys.Omega * phys.a / z0),
             0.0,
             1.0,
             duals(n) + sqr(phys.Omega * phys.a) / z0);

      hs.add_halfspace(H);
    }

#ifdef PROFILING
    time->end_section();
#endif
  }

  void extract_diagram() {
#ifdef PROFILING
    time->start_section(idx_time_extract_vertices);
#endif

    // extract vertices
    verts.resize(hs.mesh.dvert.capacity(), 2);
    dfacets.resize(hs.mesh.dvert.capacity());
    for (int i = 0; i < (int)verts.rows(); i++)
      verts.row(i) = (hs.mesh.dvert[i].head(2) / hs.mesh.dvert[i](3)).template cast<double>();

#ifdef PROFILING
    time->start_section(idx_time_extract_edges);
#endif

    // extract edges and faces
    faces.resize(hs.mesh.pcnt());
    edglist.reserve(hs.mesh.dcnt());
    for (int i = 0; i < n; i++)
      faces[i].reserve(hs.mesh.padj[i + 6].size());
    for (int i = 0; i < n; i++)
      for (const auto &e : hs.mesh.padj[i + 6])
        if (e.pj < 6 || e.pi < e.pj) {
          edglist.push_back({e.pi, e.pj, e.di, e.dj,
              discretized_line_segment(verts.row(e.di), verts.row(e.dj),
                                       phys, std::numeric_limits<double>::infinity())});

          faces[e.pi].push_back(edglist.size() - 1);
          faces[e.pj].push_back(edglist.size() - 1);
          dfacets[e.di].push_back(edglist.size() - 1);
          dfacets[e.dj].push_back(edglist.size() - 1);
        }

#ifdef PROFILING
    time->start_section(idx_time_calc_areas);
#endif

    // calculate face areas
    areas.resize(n);
    areaerrs.resize(n);
    areas.array() = 0;
    areaerrs.array() = 0;

    for (int i = 0; i < n; i++) {

      // calculate initial area and error
      for (auto &ei : faces[i + 6]) {
        auto &e = edglist[ei];
        areas(i) += (e.pi == i + 6 ? e.ls.area : -e.ls.area);
        areaerrs(i) += e.ls.errb;
      }

      if (areas(i) == 0 || areaerrs(i) / areas(i) < sim.area_tolerance) {
        // py::print(i, "- error already in tol:", areaerrs(i) / areas(i));
        continue;
      }

      std::multimap<double, int> errmap;

      for (auto &ei : faces[i + 6]) {
        auto &e = edglist[ei];
        errmap.insert({-e.ls.errb, ei});
      }

      for (int j = 0; errmap.size() > 0 && j < sim.max_refine_steps && areaerrs(i) / areas(i) >= sim.area_tolerance; j++) {
        auto [err, ei] = *errmap.begin();
        errmap.erase(errmap.begin());

        auto &e = edglist[ei];

        areas(i) -= (e.pi == i + 6 ? e.ls.area : -e.ls.area);
        areaerrs(i) -= e.ls.errb;

        e.ls.refine();

        areas(i) += (e.pi == i + 6 ? e.ls.area : -e.ls.area);
        areaerrs(i) += e.ls.errb;

        // py::print(i, "- refining edge", ei, "with errb ", e.ls.errb, "to:", areaerrs(i) / areas(i));

        errmap.insert({-e.ls.errb, ei});
      }

      // py::print(i, "- error after refining:", areaerrs(i) / areas(i));
    }

#ifdef PROFILING
    time->start_section(idx_time_calc_integrals);
#endif

    // calculate line integrals
    for (auto &e : edglist) {

      e.dif = 0;

      // points and tangents along path
      std::vector<Eigen::Vector2d> x; //, t;
      x.reserve(e.ls.points.size());
      //t.reserve(e.ls.points.size());

      for (const auto &tp : e.ls.points) {
        x.push_back(tp.x);
        // t.push_back(tp.t);
      }

      auto yi = ys.row(e.pi-6);

      if (e.pj >= 6 && e.pj - 6 < n) {
        // inner edges
        auto yj = ys.row(e.pj-6);

        double pw = 0;
        for (int k = 0; k < (int)x.size(); ++k) {
          double nw = (k != (int)x.size() - 1 ? (x[k+1] - x[k]).norm() : 0);

          // d_x (c(x,yi) - c(x,yj))
          Eigen::Vector2d dc((sqr(yi(0)) - sqr(yj(0))) * x[k](0) / (sqr(phys.a) * sqr(1 - sqr(x[k](0)))),
                             phys.cp * phys.kappa / phys.p00 * (yi(1) - yj(1)) * std::pow(x[k](1) / phys.p00, phys.kappa-1));
          e.dif += (pw + nw) / dc.norm(); // * t[k].norm() / t[k].cross(dc);

          pw = nw;
        }
      } else if (e.pj - 6 >= n) {
        // (soft) boundary edges

        double pw = 0;
        for (int k = 0; k < (int)x.size(); ++k) {
          double nw = (k != (int)x.size() - 1 ? (x[k+1] - x[k]).norm() : 0);

          // d_x (c(x,yi) - c(x,yj))
          Eigen::Vector2d dc(sqr(yi(0)) * x[k](0) / (sqr(phys.a) * sqr(1 - sqr(x[k](0)))) - sqr(phys.Omega * phys.a) * x[k](0),
                             phys.cp * phys.kappa / phys.p00 * yi(1) * std::pow(x[k](1) / phys.p00, phys.kappa-1));
          e.dif += (pw + nw) / dc.norm(); //* t[k].norm() / t[k].cross(dc);

          pw = nw;
        }

      }
    }

#ifdef PROFILING
    time->end_section();
#endif
  }

  VecX touching_dual(bool randomize) {
    // get the dual (except last coord) that touches the polytope
    // for primal coords that are part of the polytope, this is unchanged
    // for primal coords that are not part, it is larger
    // equivalent to c-transform on other dual
    VecX res = VecX::Zero(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.25, 0.75);
    
    for (int pi = 0; pi < n; ++pi) {
      if (hs.mesh.pneigh(pi + 6).size() > 0) {
        // plane already intersects -> keep current value of dual
        res[pi] = duals[pi];
      } else {
        // plane does not currently intersect...
        
        // direction of halfplane at index pi
        auto pp = hs.mesh.pvert[pi + 6].head(3);
        // z-offset of plane (without dual potential part)
        auto os = phys.Omega * ys(pi, 0);
        // find closest dual point
        int di = hs.mesh.extremal_dual(pp);
        auto dp = hs.mesh.dvert[di].head(3) / hs.mesh.dvert[di][3];
        
        auto d0 = pp.dot(dp);
        // mathematically, at phi[pi] := (-d0 - os), the plane would be touching (usually in one point)
        // however, rounding error may make it barely not touch, so we make sure it does so with a margin by looking at the neighborhood

        // get primal face that is first cut off
        auto [pj, v] = hs.mesh.closest_primal(pp);

        // set value between first touch and cutoff of primal face
        if (randomize) {
          double w = dis(gen);
          res[pi] = -(w * d0 + (1 - w) * v) - os;
        } else {
          res[pi] = -0.5 * (d0 + v) - os;
        }
      }
    }

    return res;
  }

  rasterizer get_rasterizer(const double &merge_epsi) {

#ifdef PROFILING
    time->start_section(idx_time_prepare_rasterizer);
#endif

    // prepare segment indices
    std::vector<int> sidx(edglist.size() + 1); // index of first segment for each edge
    sidx[0] = 0;
    for (int i = 0; i < (int)edglist.size(); ++i) {
      const auto &e = edglist[i];

      // primal indices 1 == bottom, 2 == left, 3 == right, 4 == top
      if (e.pj == 2 || e.pj == 3) // skip vertical segments on left and right
        sidx[i + 1] = sidx[i];
      else
        sidx[i + 1] = sidx[i] + e.ls.points.size() - 1;
    }

    // rasterizer arguments
    std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> seg;
    std::vector<std::vector<int>> con;
    std::vector<int> start;
    Eigen::Array4d bounds = {sim.spmin(0), sim.spmax(0), sim.spmin(1), sim.spmax(1)};

    // collect segments and transition connections
    for (int i = 0; i < (int)edglist.size(); ++i) {
      const auto &e = edglist[i];
      const auto &points = e.ls.points;

      if (e.pj == 2 || e.pj == 3)
        continue;

      // get (primal) index of cell below:
      // if the edge is oriented left-to-right then pi is below, if it goes right-to-left then pj is below
      // indices >= n (i.e. the top cell) are all combined into n
      int vidx = std::min(((points.front().x(0) > points.back().x(0)) ? e.pj : e.pi) - 6, n);
      if (vidx < 0) vidx = std::min(((points.front().x(0) > points.back().x(0)) ? e.pi : e.pj) - 6, n);
      if (vidx < 0) throw std::runtime_error("in preparing rasterizer, edge has both primal vertices outside: " + e.repr());

      // add segments
      int j = sidx[i];
      for (auto it1 = points.begin(), it0 = it1++; it1 != points.end(); it0 = it1++, ++j) {
        if ((int)seg.size() != j)
          throw std::runtime_error(FORMAT("segment index does not match expected: {} instead of {}", seg.size(), j));
        seg.push_back({it0->x, it1->x, vidx});
      }

      // add transition connections
      for (int j = 1; j < (int)points.size() - 1; ++j)
        con.push_back({sidx[i] + j - 1, sidx[i] + j});
    }

    // collect merge/split points
    for (int dc = 0; dc < (int)dfacets.size(); ++dc) {
      const auto &f = dfacets[dc];

      // ignore empty facets
      if (!f.size())
        continue;

      // find simulation box boundaries the point is on
      bool on_boundary[6] = {};
      for (const int &ei : f)
        if (edglist[ei].pj < 6)
          on_boundary[edglist[ei].pj] = true;

      // check for top / bottom connections (should not exist)
      if (on_boundary[0] || on_boundary[5])
        throw std::runtime_error("found dual point touching 3d limit");

      // if the point is on the left boundary, add the edge going away to start
      if (on_boundary[2]) { // x0 lo (left)
        for (const int &ei : f)
          if (edglist[ei].pj != 2)
            start.push_back(edglist[ei].di == dc ? sidx[ei] : sidx[ei + 1] - 1);
      }

      // if the point is not on the left or right boundary, add transition / merge / split connection
      if (!on_boundary[2] && !on_boundary[3]) {
        std::vector<int> c;
        for (const int &ei : f)
          c.push_back(edglist[ei].di == dc ? sidx[ei] : sidx[ei + 1] - 1);
        con.push_back(c);
      }
    }

    seg.push_back({{sim.spmin(0), sim.spmax(1) + 1}, {sim.spmax(0), sim.spmax(1) + 1}, n});
    start.push_back(seg.size() - 1);

    // seg.push_back({{sim.spmin(0), sim.spmin(1) - 1}, {sim.spmax(0), sim.spmin(1) - 1}, n});
    // start.push_back(seg.size() - 1);

    rasterizer rast(seg, con, start, bounds, merge_epsi);

#ifdef PROFILING
    time->end_section();
#endif

    return rast;
  }

  std::unordered_map<uint64_t, double> jac() {
#ifdef PROFILING
    time->start_section(idx_time_jacobian);
#endif
    // get jacobian (derivative of masses w.r.t. dual)

    std::unordered_map<uint64_t, double> res;
    for (auto &e : edglist) {
      if (e.pj >= 6) {
        uint64_t i = e.pi-6;
        uint64_t j = std::min(e.pj-6,n);

        res[(i << 32) | i] += e.dif;
        res[(i << 32) | j] -= e.dif;

        res[(j << 32) | j] += e.dif;
        res[(j << 32) | i] -= e.dif;
      }
    }
#ifdef PROFILING
    time->end_section();
#endif
    return res;
  }

  std::pair<Eigen::VectorXd, std::pair<Eigen::VectorXi, Eigen::VectorXi>> jac_coo() {
    auto J = jac();

#ifdef PROFILING
    time->start_section(idx_time_jacobian);
#endif

    Eigen::VectorXd data(J.size());
    Eigen::VectorXi i(J.size()), j(J.size());

    int k = 0;
    for (auto &v : J) {
      data[k] = v.second;
      i[k] = v.first >> 32;
      j[k++] = v.first & 0xFFFFFFFF;
      // i[k] = v.first.first;
      // j[k++] = v.first.second;
    }

#ifdef PROFILING
    time->end_section();
#endif
    return {data, {i, j}};
  }

  Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> get_poly(const int &pi_t) {
    int pi = pi_t + 6;
    std::vector<diagram_edge*> edg, path;
    int len = 0;
    for (const auto &ei : faces[pi]) {
      edg.push_back(&edglist[ei]);
      len += edglist[ei].ls.points.size() - 1;
    }

    if (edg.size() <= 2)
      throw std::runtime_error("polygon only has 2 sides");

    std::map<int, std::vector<int>> adj;
    for (int i = 0; i < (int)edg.size(); ++i) {
      adj[edg[i]->pi == pi ? edg[i]->di : edg[i]->dj].push_back(i);
    }
    std::vector<bool> used(edg.size(), 0);

    std::function<void(const int&)> dfs = [&](const int &ev){
      int v = edg[ev]->pi == pi ? edg[ev]->dj : edg[ev]->di;
      for (const int &ew : adj[v])
        if (!used[ew])
          used[ew] = true, dfs(ew);
      path.push_back(edg[ev]);
    };
    used[0] = 1;
    dfs(0);
    std::reverse(path.begin(), path.end());

    if (path.size() != edg.size())
      throw std::runtime_error(FORMAT("failed to orient polygon {}: got {} edges of {}", pi, path.size(), edg.size()));
    
    int ri = 0;
    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> res(len, 2);
    for (auto &ep : path)
      if (ep->pi == pi)
        for (auto pit = ep->ls.points.begin(); pit != std::prev(ep->ls.points.end()); ++pit)
          res.row(ri++) = pit->x;
      else
        for (auto pit = std::prev(ep->ls.points.end()); pit != ep->ls.points.begin(); --pit)
          res.row(ri++) = pit->x;

    return res;
  }
};

#define BIND_DIAGRAM_EDGE(m, T, name)                                   \
  py::class_<laguerre_diagram<T>::diagram_edge>(m, name)                \
  .def("__repr__", &laguerre_diagram<T>::diagram_edge::repr)            \
  .def_readonly("pi", &laguerre_diagram<T>::diagram_edge::pi)           \
  .def_readonly("pj", &laguerre_diagram<T>::diagram_edge::pj)           \
  .def_readonly("di", &laguerre_diagram<T>::diagram_edge::di)           \
  .def_readonly("dj", &laguerre_diagram<T>::diagram_edge::dj)           \
  .def_readonly("ls", &laguerre_diagram<T>::diagram_edge::ls)           \
  .def_readonly("dif", &laguerre_diagram<T>::diagram_edge::dif);


#ifdef PROFILING
#define BIND_LAGUERRE_DIAGRAM_PROFILING(m, T)           \
  def_readonly("time", &laguerre_diagram<T>::time)
#else
#define BIND_LAGUERRE_DIAGRAM_PROFILING(m, T)
#endif

#define BIND_LAGUERRE_DIAGRAM(m, T, name)                       \
  py::class_<laguerre_diagram<T>>(m, name)                      \
  .def(py::init<const laguerre_diagram<T>::seeds_t&,            \
       const Eigen::Ref<const laguerre_diagram<T>::VecX>&,      \
       const physical_parameters &,                             \
       const simulation_parameters&>(),                         \
       py::arg("ys"), py::arg("duals"),                         \
       py::arg("phys"), py::arg("sim"))                         \
  .def_readonly("n", &laguerre_diagram<T>::n)                   \
  .def_readonly("phys", &laguerre_diagram<T>::phys)             \
  .def_readonly("sim", &laguerre_diagram<T>::sim)               \
  .def_readonly("ys", &laguerre_diagram<T>::ys)                 \
  .def_readonly("duals", &laguerre_diagram<T>::duals)           \
  .def_readonly("hs", &laguerre_diagram<T>::hs)                 \
  .def_readonly("edglist", &laguerre_diagram<T>::edglist)       \
  .def_readonly("verts", &laguerre_diagram<T>::verts)           \
  .def_readonly("faces", &laguerre_diagram<T>::faces)           \
  .def_readonly("dfacets", &laguerre_diagram<T>::dfacets)       \
  .def_readonly("areas", &laguerre_diagram<T>::areas)           \
  .def_readonly("areaerrs", &laguerre_diagram<T>::areaerrs)     \
  .BIND_LAGUERRE_DIAGRAM_PROFILING(m, T)                        \
  .def("touching_dual", &laguerre_diagram<T>::touching_dual,    \
       py::arg("randomize")=false)                              \
  .def("jac_coo", &laguerre_diagram<T>::jac_coo)                \
  .def("get_rasterizer", &laguerre_diagram<T>::get_rasterizer,  \
       py::arg("merge_epsi") = 0)                               \
  .def("get_poly", &laguerre_diagram<T>::get_poly);


#endif
