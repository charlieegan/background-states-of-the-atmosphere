#ifndef HALFSPACE_INTERSECTION_HPP
#define HALFSPACE_INTERSECTION_HPP

#include "common.hpp"
#include "dvector.hpp"
#include "reusable_used_array.hpp"
#include "timer.hpp"

/*! An (oriented) primal-dual edge represented by primal and dual indices. */
struct pdedge {
  int pi; //!< primal starting index
  int pj; //!< primal ending index
  int di; //!< dual starting index
  int dj; //!< dual ending index

  /*! Return a string representation of *this. */
  std::string repr() const {
    std::stringstream ss;
    ss << "PDEdge(" << pi << ", " << pj << "; " << di << ", " << dj << ")";
    return ss.str();
  }

  /*! Get mirrored edge (same orientation but opposite direction). */
  pdedge mirrored() const {
    return {pj, pi, dj, di};
  }

  /*! Functor struct to compare pdedge by primal index. */
  struct cmp_pj { bool operator()(const pdedge &e, const pdedge &f) const { return e.pj < f.pj; } };

  /*! Functor struct to compare pdedge by dual index. */
  struct cmp_dj { bool operator()(const pdedge &e, const pdedge &f) const { return e.dj < f.dj; } };

  /*! Order the given edges to make up a cycle on primal side (if primal == true) or dual side (if primal == false).
   *  If the edges do not make up a cycle, a runtime_error will the thrown, the function assumes that edges are already
   *  oriented the same way and form a cycle under a suitable permutation.
   * \param edg vector of input edges
   * \return permutation of edg such that r[i].pj == r[(i+1)%n].pi (if primal == true, dual correspondingly)
   */
  template <bool primal = true>
  static std::vector<pdedge> order_cycle(const std::vector<pdedge> &edg) {
    if (edg.size() <= 2)
      return edg;

    std::map<int, std::vector<int>> adj;
    for (int i = 0; i < (int)edg.size(); ++i) {
      adj[primal ? edg[i].pi : edg[i].di].push_back(i);
    }
    std::vector<bool> used(edg.size(), 0);
    std::vector<pdedge> path;

    std::function<void(const int&)> dfs = [&](const int &ev){
      int v = (primal ? edg[ev].pj : edg[ev].dj);

      for (const int &ew : adj[v]) {
        if (!used[ew]) {
          // int w = (primal ? edg[ew].pj : edg[ew].dj);

          used[ew] = true;
          dfs(ew);
        }

      }
      path.push_back(edg[ev]);
    };

    used[0] = 1;
    dfs(0);

    std::reverse(path.begin(), path.end());

#ifdef DEBUG_CHECKS
    if (path.size() != edg.size())
      throw std::runtime_error("resulting path length in order_cycle is too short");

    if (primal ? (path.back().pj != path.front().pi) : (path.back().dj != path.front().di))
      throw std::runtime_error("path in order_cycle does not close");
#endif

    return path;
  }

  /*! Add python bindings of pdedge to m. */
  static void bind(py::module_ &m) {
    py::class_<pdedge>(m, "PDEdge")
      .def("__repr__", &pdedge::repr)
      .def_readonly("pi", &pdedge::pi)
      .def_readonly("pj", &pdedge::pj)
      .def_readonly("di", &pdedge::di)
      .def_readonly("dj", &pdedge::dj);
  }
};

/*! An ordered primal-dual mesh in 3d (projective) space with dual vertices having degree 3. */
template <typename T = double>
struct pdmesh {
  typedef Eigen::Vector3<T> Vec3;
  typedef Eigen::Vector4<T> Vec4;

  std::vector<std::vector<pdedge>> padj; //!< primal adjacency list - for each primal index: list of incident edges
  dvector<std::vector<pdedge>> dadj; //!< dual adjacency list - for each dual index: list of incident edges (length 3)

  std::vector<Vec4> pvert; //!< primal vertices
  dvector<Vec4> dvert; //!< dual vertices

  reusable_used_array<> used; //!< used array for internal algorithms

  /*! Construct mesh allocating space for primal and dual vertices if requested.
   * \param pcnt number of primal vertices to allocate space for
   * \param dcnt number of dual vertices to allocate space for
   */
  pdmesh(const int &pcnt = 0, const int &dcnt = 0) :
    padj(pcnt), dadj(dcnt),
    pvert(pcnt), dvert(dcnt) {}

  /*! Return a cuboid pdmesh with given lower and upper bounds.
   * \param lb corner with lowest coordinates
   * \param ub corner with highest coordinates
   * \return pdmesh representing the cuboid with given dimensions
   */
  static pdmesh cube(const Eigen::Ref<const Vec3> &lb,
                     const Eigen::Ref<const Vec3> &ub) {
    pdmesh res(6, 8);

    // a cuboid mesh has 6 primal vertices, 8 dual vertices, 12 edges

    /* enumeration of primal and dual vertices is as follows:
     *
     *      4-----------5
     *     /|          /|
     *    / |   5     / |       x2
     *   /  |     1  /  |       |
     *  6-----------7   |       |
     *  | 2 |       | 3 |       |
     *  |   0-------|---1       o ----- x0
     *  |  /  4     |  /       /
     *  | /     0   | /       /
     *  |/          |/       x1
     *  2-----------3
     */

    // define primal vertex (i.e. halfspace) coordinates
    res.pvert[0] = { 0,  0, -1,  lb(2)};
    res.pvert[5] = { 0,  0,  1, -ub(2)};
    res.pvert[1] = { 0, -1,  0,  lb(1)};
    res.pvert[4] = { 0,  1,  0, -ub(1)};
    res.pvert[2] = {-1,  0,  0,  lb(0)};
    res.pvert[3] = { 1,  0,  0, -ub(0)};

    res.dvert[0] = {lb[0], lb[1], lb[2], 1.0};
    res.dvert[1] = {ub[0], lb[1], lb[2], 1.0};
    res.dvert[2] = {lb[0], ub[1], lb[2], 1.0};
    res.dvert[3] = {ub[0], ub[1], lb[2], 1.0};
    res.dvert[4] = {lb[0], lb[1], ub[2], 1.0};
    res.dvert[5] = {ub[0], lb[1], ub[2], 1.0};
    res.dvert[6] = {lb[0], ub[1], ub[2], 1.0};
    res.dvert[7] = {ub[0], ub[1], ub[2], 1.0};

    // sign convention:
    // when standing on the inside of the domain, looking from pi over the dual edge to pj
    // then di is on the left and dj is on the right.

    // when following an the edges starting at a dual vertex, e.g. for di = 0
    // {0, 1, _0_, 1}, {1, 2, _0_, 4}, {2, 0, _0_, 2}
    // looking from the inside, they go around counter-clockwise

    // bottom edges
    res.add_edge({0, 1, 0, 1});
    res.add_edge({0, 2, 2, 0});
    res.add_edge({0, 3, 1, 3});
    res.add_edge({0, 4, 3, 2});

    // side edges
    res.add_edge({2, 1, 4, 0});
    res.add_edge({1, 3, 5, 1});
    res.add_edge({3, 4, 7, 3});
    res.add_edge({4, 2, 6, 2});

    // top edges
    res.add_edge({5, 2, 4, 6});
    res.add_edge({5, 4, 6, 7});
    res.add_edge({5, 3, 7, 5});
    res.add_edge({5, 1, 5, 4});

#ifdef DEBUG_CHECKS
    if (!res.check_integrity())
      throw std::runtime_error("generated cube has wrong integrity");
#endif

    return res;
  }

  /*! Get number of primal indices. */
  int pcnt() const {
#ifdef DEBUG_CHECKS
    if (padj.size() != pvert.size())
      throw std::runtime_error("padj and pvert do not have the same size");
#endif
    return padj.size();
  }

  /*! Get number of dual indices. */
  int dcnt() const {
#ifdef DEBUG_CHECKS
    if (dadj.size() != dvert.size())
      throw std::runtime_error("dadj and dvert do not have the same size");
#endif
    return dadj.size();
  }

  /*! Add primal-dual edge to mesh. */
  void add_edge(const pdedge &e) {
#ifdef DEBUG_CHECKS
    if (e.pi < 0 || e.pj < 0 || e.pi >= pcnt() || e.pj >= pcnt() ||
        !dadj.contains_idx(e.di) || !dadj.contains_idx(e.dj))
      throw std::runtime_error("tried to add edge containing invalid points " + e.repr());
#endif

    padj[e.pi].push_back(e);
    padj[e.pj].push_back(e.mirrored());
    dadj[e.di].push_back(e);
    dadj[e.dj].push_back(e.mirrored());

    fix_order(e.di);
    fix_order(e.dj);
  }

  /*! Remove primal-dual edge from mesh. */
  void remove_edge(const pdedge &e) {
    // TODO: efficiency (using vectors may actually be best when number of elements is very small -> benchmark)
#ifdef DEBUG_CHECKS
    int i0 =
#endif
      std::erase_if(padj[e.pi], [&e](const pdedge &f){ return f.pj == e.pj && f.di == e.di && f.dj == e.dj; });
#ifdef DEBUG_CHECKS
    int i1 =
#endif
      std::erase_if(padj[e.pj], [&e](const pdedge &f){ return f.pj == e.pi && f.di == e.dj && f.dj == e.di; });
#ifdef DEBUG_CHECKS
    int i2 =
#endif
      std::erase_if(dadj[e.di], [&e](const pdedge &f){ return f.dj == e.dj && f.pi == e.pi && f.pj == e.pj; });
#ifdef DEBUG_CHECKS
    int i3 =
#endif
      std::erase_if(dadj[e.dj], [&e](const pdedge &f){ return f.dj == e.di && f.pi == e.pj && f.pj == e.pi; });

#ifdef DEBUG_CHECKS
    if (i0 < 1 || i1 < 1 || i2 < 1 || i3 < 1)
      throw std::runtime_error("tried to erase edge that does not exist");
    if (i0 > 1 || i1 > 1 || i2 > 1 || i3 > 1)
      throw std::runtime_error("erase removed multiple edges");
#endif

    fix_order(e.di);
    fix_order(e.dj);
  }

  /*! Remove dual vertex at index with all corresponding edges. */
  void remove_dual(const int &di) {
#ifdef DEBUG_CHECKS
    if (!dadj.contains_idx(di))
      throw std::runtime_error("tried to erase dual that does not exist");
#endif

    while (dadj[di].size()) {
      remove_edge(dadj[di].back());
    }
    dadj.remove(di);
    dvert.remove(di);
  }

  /*! Add a new dual, return the index. */
  int add_dual(const Eigen::Ref<const Vec4> &p) {
    dvert.add(p);
    int di = dadj.add(std::vector<pdedge>());
    return di;
  }

  /*! Add a new primal vertex and return the index. */
  int add_primal(const Eigen::Ref<const Vec4> &hs) {
    pvert.push_back(hs);
    padj.push_back(std::vector<pdedge>());
    return (int)pvert.size() - 1;
  }

  /*! Make sure the order of dadj[di] is correct if dadj[di].size() == 3. */
  void fix_order(const int &di) {
#ifdef DEBUG_CHECKS
    if (!dadj.contains_idx(di))
      throw std::runtime_error("fix_order called with invalid index");
#endif

    if (dadj[di].size() != 3)
      return;

    if (dadj[di][0].pj != dadj[di][1].pi)
      std::swap(dadj[di][1], dadj[di][2]);

#ifdef DEBUG_CHECKS
    if (dadj[di][0].pj != dadj[di][1].pi ||
        dadj[di][1].pj != dadj[di][2].pi ||
        dadj[di][2].pj != dadj[di][0].pi)
      throw std::runtime_error("dual edges adjacent in dadj do not form cycle");
#endif

    // auto &a = pvert[dadj[di][0].pi];
    // auto &b = pvert[dadj[di][1].pi];
    // auto &c = pvert[dadj[di][2].pi];

    // dvert[di] =
    //   {
    //     -a({1,2,3}).cross(b({1,2,3})).dot(c({1,2,3})),
    //     a({2,3,0}).cross(b({2,3,0})).dot(c({2,3,0})),
    //     -a({3,0,1}).cross(b({3,0,1})).dot(c({3,0,1})),
    //     a({0,1,2}).cross(b({0,1,2})).dot(c({0,1,2}))
    //   };

#ifdef DEBUG_CHECKS
    if (dvert[di](3) <= 1e-9) {
      throw std::runtime_error(FORMAT("halfspaces in dual intersection vertex {}, intersect of primals {}, {}, {} not geometrically ordered: det = {}",
                                      di, dadj[di][0].pi, dadj[di][1].pi, dadj[di][2].pi, dvert[di](3)));
    }
#endif

  }

  /*! Find how far inside / outside of the halfspace corresponding to pi the dual point di is.
   *  The value is negative if the dual point is inside the halfspace, positive if it is outside.
   */
  T hs_error(const int &pi, const int &di) const {
    return pvert[pi].dot(dvert[di]) / (dvert[di](3) * pvert[pi].head(3).norm());
  }

  /*! Check whether halfspace defined by pi contains dual point di. */
  bool contains(const int &pi, const int &di) const {
#ifdef DEBUG_CHECKS
    if (pi < 0 || pi >= pcnt() || !dadj.contains_idx(di))
      throw std::runtime_error("called contains() with invalid points");
#endif
    return pvert[pi].dot(dvert[di]) <= 0;
  }

  /*! Find the index of the dual point which has largest inner product with d. */
  int extremal_dual(const Eigen::Ref<const Vec3> &d, bool brute_force=false, int hint = -1, bool debug=false) {
    int res = (hint >= 0 ? hint : *dadj.begin_idx());
    T maxval = dvert[res].head(3).dot(d) / dvert[res](3), val;

    if (brute_force) {
      for (auto it = dadj.begin_idx(); it != dadj.end_idx(); ++it) {
        int di = *it;
        val = dvert[di].head(3).dot(d) / dvert[di](3);
        if (val > maxval) {
          res = di;
          maxval = val;
        }
      }
    } else {
      int tcnt = 0;

      // TODO: optimze this value or at least make it accessible
      const T eps = 1e-5;

      // traverse mesh to find extremal point
      used.reset();
      std::deque<std::pair<T, int>> q;
      used.mark(res);
      q.push_back({maxval, res});
      while (!q.empty()) {
        auto [v, di] = q.back(); q.pop_back();
        tcnt++;

        // go over neighbors and push them if they are not too far back and not already visited
        for (const auto &e : dneigh(di)) {
          val = dvert[e.dj].head(3).dot(d) / dvert[e.dj](3);
          if (!significantly_less<T>(val, maxval, eps, eps) && !used.check(e.dj)) {

            used.mark(e.dj);
            if (val > maxval)
              q.push_back({val, e.dj});
            else
              q.push_front({val, e.dj});

            // if we found new max, store it
            if (val > maxval) {
              maxval = val;
              res = e.dj;
            }
          }
        }
      }
      if (debug)
        py::print("tcnt:", tcnt, "pqlen:", q.size());
    }

#ifdef DEBUG_CHECKS
    for (const auto &e : dneigh(res)) {
      val = dvert[e.dj].head(3).dot(d) / dvert[e.dj](3);
      if (val > maxval)
        throw std::runtime_error("extremal dual is not locally extremal");
    }
#endif

    return res;
  }

  /*! Find the (index of the) primal point which is first cut off by the plane with normal d when translating it.
   * \return {pi, <d, last (dual) vertex of pi cut off>}
   */
  std::pair<int, T> closest_primal(const Eigen::Ref<const Vec3> &d) {

    T maxval = -std::numeric_limits<T>::infinity(), val;
    int res = -1;

    for (int pi = 6; pi < (int)padj.size(); ++pi) {

      if (padj[pi].size() == 0)
        continue;

      val = std::numeric_limits<T>::infinity();
      for (const auto &e : padj[pi]) {
        val = std::min(val, dvert[e.di].head(3).dot(d) / dvert[e.di](3));
        if (val < maxval)
          break;
      }

      if (val > maxval) {
        res = pi;
        maxval = val;
      }
    }

    return {res, maxval};
  }

  /*! Fet neighbors of primal index pi (as edges e with e.pi == pi). */
  const std::vector<pdedge>& pneigh(const int &pi) const {
#ifdef DEBUG_CHECKS
    if (pi < 0 || pi >= pcnt())
      throw std::runtime_error("called pneigh with invalid index");
#endif
    return padj[pi];
  }

  /*! Get neighbors of dual index di (as edges e with e.di == di). */
  const std::vector<pdedge>& dneigh(const int &di) const {
#ifdef DEBUG_CHECKS
    if (!dadj.contains_idx(di))
      throw std::runtime_error("called pneigh with invalid index");
#endif
    return dadj[di];
  }

  /*! Get a minimalistic string representation. */
  std::string repr() const {
    std::stringstream ss;
    ss << "PDMesh(pcnt = " << pcnt() << ", dcnt = " << dcnt() << ")";
    return ss.str();
  }

  /*! Get the polygon (list of dual vertices) corresponding to primal index pi. */
  std::vector<Vec4> primal_poly(const int &pi) const {
    std::vector<Vec4> res;
    for (const auto &e : pdedge::order_cycle<false>(padj[pi])) {
      res.push_back(dvert[e.di]);
    }
    return res;
  }

  /*! Do a bunch of checks for mesh integrity (for debugging purposes). */
  bool check_integrity() const {
    bool res = true;

    if (dadj.size() != dvert.size()) {
      py::print("dadj and dvert have different sizes");
      res = false;
    }
    if (padj.size() != pvert.size()) {
      py::print("padj and pvert have different sizes");
      res = false;
    }

    for (const auto &v : dadj) {
      // check dual face size
      if (v.size() != 3) {
        py::print("non-triangular dual face found:", v);
        res = false;
      }

      // check dual face orientation / order
      int first_pi = -1, prev_pj = -1;
      bool first = true;
      for (const auto &e : v) {
        if (first) {
          first_pi = e.pi;
          first = false;
        } else {
          if (e.pi != prev_pj) {
            py::print("incorrectly ordered or oriented dual face found:", v);
            res = false;
          }
        }
        prev_pj = e.pj;
      }
      if (first_pi != prev_pj) {
        py::print("incorrectly ordered dual face found: ", v);
        res = false;
      }
    }
    for (const auto &v : padj) {
      // check primal face size
      if (v.size() == 1 || v.size() == 2) {
        py::print("degenerated primal face found:", v);
        res = false;
      }
    }

    for (int pi = 0; pi < (int)pvert.size(); ++pi)
      for (auto dit = dvert.begin_idx(); dit != dvert.end_idx(); ++dit) {
        int di = *dit;
        auto v = pvert[pi];
        auto dv = dvert[di];

        bool neigh = false;
        for (const auto &e : padj[pi])
          if (di == e.di)
            neigh = true;

        T dd = v.dot(dv);
        if (!neigh && dd > 1e-9) {
          py::print(FORMAT("convexity violation by {} @ primal {} ({}, {}, {}, {}), dual {} ({}, {}, {}, {})",
                           dd, pi, v(0), v(1), v(2), v(3), di, dv(0), dv(1), dv(2), dv(3)));
          res = false;
        }
      }

    return res;
  }

  /*! Write mesh to a .ply file at given path (for debugging purposes). */
  void write_ply(std::string path) const {
    std::ofstream s(path);
    s << std::setprecision(10);
    s << "ply" << std::endl
      << "format ascii 1.0" << std::endl
      << "element vertex " << dvert.capacity() << std::endl
      << "property float x" << std::endl
      << "property float y" << std::endl
      << "property float z" << std::endl
      << "element face " << pcnt() << std::endl
      << "property list uint uint vertex_indices" << std::endl
      << "end_header" << std::endl;
    for (const auto &v : dvert.data) {
      s << std::clamp(v[0] / v[3], (T)-1e30, (T)1e30) << " "
        << std::clamp(v[1] / v[3], (T)-1e30, (T)1e30) << " "
        << std::clamp(v[2] / v[3], (T)-1e30, (T)1e30) << std::endl;
    }
    for (const auto &es : padj) {
      s << es.size();
      try {
        for (const auto &e : pdedge::order_cycle<false>(es))
          s << " " << e.di;
      } catch (std::runtime_error &ex) {
        py::print("WARNING: failed to order_cycle in write_ply");
        for (const auto &e : es)
          s << " " << e.di;
      }
      s << std::endl;
    }
  }

  /*! Fix order of all primal faces. */
  void order_all() {
    for (int i = 0; i < (int)padj.size(); ++i) {
      padj[i] = pdedge::order_cycle<false>(padj[i]);
    }
  }

  /*! Add python binding code for module m. */
  static void bind(py::module_ &m) {
    py::class_<pdmesh<T>>(m, ("PDMesh_" + type_name<T>::value()).c_str())
      .def("__repr__", &pdmesh<T>::repr)
      .def("check_integrity", &pdmesh<T>::check_integrity)
      .def("primal_poly", &pdmesh<T>::primal_poly)
      .def_property_readonly("pcnt", &pdmesh<T>::pcnt)
      .def_property_readonly("dcnt", &pdmesh<T>::dcnt)
      .def_readonly("padj", &pdmesh<T>::padj)
      .def_property_readonly("dadj",
                             [](const pdmesh<T> *m){ return m->dadj.data;})
      .def_readonly("pvert", &pdmesh<T>::pvert)
      .def_property_readonly("dvert",
                             [](const pdmesh<T> *m){ return m->dvert.data;})
      .def("extremal_dual", &pdmesh<T>::extremal_dual,
           py::arg("d"), py::arg("brute_force")=false,
           py::arg("hint")=-1, py::arg("debug")=false)
      .def("closest_primal", &pdmesh<T>::closest_primal)
      .def("write_ply", &pdmesh<T>::write_ply)
      .def("order_all", &pdmesh<T>::order_all);
  }
};

template <typename T = double>
class halfspace_intersection
{
public:
  typedef Eigen::Vector3<T> Vec3;
  typedef Eigen::Vector4<T> Vec4;

  pdmesh<T> mesh; //!< primal-dual mesh representing intersection of halfspaces so far
  reusable_used_array<> used; //!< used array for internal graph algorithms

#ifdef PROFILING
  std::shared_ptr<timer> time = NULL; //!< profiling timer
  int idx_time_find_cut = -1, idx_time_partition = -1, idx_time_cycledec = -1, idx_time_add = -1, idx_time_remove = -1;
  int remove_count = 0;

  void setup_timer() {
    if (!time)
      time = std::make_shared<timer>();
    idx_time_find_cut = time->get_index_from_name("halfspace_intersection.find_cut");
    idx_time_partition = time->get_index_from_name("halfspace_intersection.partition");
    idx_time_cycledec = time->get_index_from_name("halfspace_intersection.cycledec");
    idx_time_add = time->get_index_from_name("halfspace_intersection.add");
    idx_time_remove = time->get_index_from_name("halfspace_intersection.remove");
  }
#endif

  /*! Construct halfspace intersection with an initial cuboid mesh with given corners. */
  halfspace_intersection(const Eigen::Ref<const Vec3> &lb = Vec3(0,0,0),
                         const Eigen::Ref<const Vec3> &ub = Vec3(1,1,1)) :
    mesh(pdmesh<T>::cube(lb, ub)), used(8) {
#ifdef DEBUG_CHECKS
    if ((lb.array() >= ub.array()).any())
      throw std::runtime_error("lb !<= ub");
#endif
  }

  /*! Construct halfspace intersection with initial cuboid and halfspaces. */
  halfspace_intersection(const Eigen::Ref<const Vec3> &lb,
                         const Eigen::Ref<const Vec3> &ub,
                         const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 4, Eigen::RowMajor>> &hss) :
    halfspace_intersection(lb, ub) {
    add_halfspaces(hss);
  }

  /*! Add list of halfspaces to intersection. */
  void add_halfspaces(const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 4, Eigen::RowMajor>> &hss) {
    for (const auto &hs : hss.rowwise())
      add_halfspace(hs);
  }

  /*! Add a halfspace to the intersection.
   * \param hs vector describing halfspace: {x | hs.head(3).dot(x) + hs[3] <= 0}
   * \param hint primal index close to which the new face will likely intersect (to speed up search)
   * \return hint for next execution - if halfspace_intersection is constructed in the same way multiple times with similar inputs, supplying the hint returned from the previous time will likely speed up the current one.
   */
  int add_halfspace(const Eigen::Ref<const Vec4> &hs, int hint = -1) {
#ifdef EXPENSIVE_DEBUG_CHECKS
    if (!mesh.check_integrity())
      throw std::runtime_error("mesh integrity broken before adding new halfspace");
#endif

#ifdef PROFILING
    if (idx_time_find_cut < 0)
      setup_timer();

    time->start_section(idx_time_find_cut);
#endif

#ifdef DEBUG_CHECKS
    if (hs.array().isNaN().any())
      throw std::runtime_error("halfspace to add contains NaN");
#endif

    // find best dual hint from primal hint
    int dhint = -1;
    if (hint >= 0 && hint < (int)mesh.padj.size() && mesh.pneigh(hint).size() > 0) {
      T bval = -std::numeric_limits<T>::infinity();
      for (const auto &e : mesh.pneigh(hint)) {
        T val = mesh.dvert[e.di].head(3).dot(hs.head(3)) / mesh.dvert[e.di][3];
        if (val > bval) {
          bval = val;
          dhint = e.di;
        }
      }
    }

#ifdef DEBUG_CHECKS
    if (dhint >= 0 && !mesh.dadj.contains_idx(dhint))
      throw std::runtime_error("dual hint is invalid");
#endif
    
    // find dual vertex furthest in direction of halfspace normal
    int di_start = mesh.extremal_dual(hs.head(3), false, dhint);

#ifdef DEBUG_CHECKS
    if (!mesh.dadj.contains_idx(di_start))
      throw std::runtime_error("extremal_dual returned invalid index");
#endif

    // find best hint for next time (since there are only 3 dual neighbors but there may be many primal, it makes sense to look for the one with smallest degree)
    int phint = -1;
    {
      int best = 1e9;
      for (const auto &e : mesh.dneigh(di_start)) {
        if ((int)mesh.pneigh(e.pi).size() < best) {
          phint = e.pi;
          best = mesh.pneigh(e.pi).size();
        }
      }
    }

#ifdef EXPENSIVE_DEBUG_CHECKS
    int di_start_brute = mesh.extremal_dual(hs.head(3), true);
    T v_trav = mesh.dvert[di_start].head(3).dot(hs.head(3)) / mesh.dvert[di_start][3];
    T v_brute = mesh.dvert[di_start_brute].head(3).dot(hs.head(3)) / mesh.dvert[di_start_brute][3];

    if (v_trav != v_brute)
      throw std::runtime_error(FORMAT("extremal_dual incorrect - trav: {} ({}, {}, {}, {}) dist ({}), "
                                      "brute: {} ({}, {}, {}, {}) dist ({}), integrity: {}",
                                      di_start, mesh.dvert[di_start](0), mesh.dvert[di_start](1), mesh.dvert[di_start](2), mesh.dvert[di_start](3), v_trav,
                                      di_start_brute, mesh.dvert[di_start_brute](0), mesh.dvert[di_start_brute](1), mesh.dvert[di_start_brute](2), mesh.dvert[di_start_brute](3), v_brute,
                                      mesh.check_integrity()));
#endif

    // add new primal
    int pi = mesh.add_primal(hs);

    // if new halfspace does not intersect, do not process further
    if (mesh.contains(pi, di_start))
      di_start = -1;

#ifdef PROFILING
    time->end_section();
#endif

    // if there is no point outside, the new halfspace does nothing
    if (di_start < 0)
      return phint;

#ifdef PROFILING
    time->start_section(idx_time_partition);
#endif

    used.reset();
    used.mark(di_start);
    std::vector<int> remove_di({di_start});  // dual vertices that are cut off (contains(pi, *) is false)
    std::vector<pdedge> cross_e;             // edges crossing from cut off part to remaining part (e.di is removed, e.dj is remaining)
    std::map<int, int> start_e;              // map from primal index pi to index in cross_e of edge with e.pi == pi

    // traverse set of cut off vertices and store the vertices and edges connecting them to the remaining part
    std::vector<int> q({di_start});
    while (!q.empty()) {
      int di = q.back();
      q.pop_back();

      for (const auto &e : mesh.dneigh(di)) {
        // ignore already visited vertices
        if (used.check(e.dj))
          continue;

        if (!mesh.contains(pi, e.dj)) {
          // e.dj is also cut off -> store and push onto queue
          q.push_back(e.dj);
          remove_di.push_back(e.dj);
          // mark used
          used.mark(e.dj);
        } else {
          // e.dj is not cut off -> store the edge
#ifdef DEBUG_CHECKS
          if (start_e.count(e.pi))
            throw std::runtime_error("multiple edges with same pi in cross_e");
#endif
          start_e[e.pi] = cross_e.size();
          cross_e.push_back(e);
        }
      }
    }

    // - remove cut dual vertices with all their edges
    // - add new dual vertices at intersections
    // - add new edges where cross_e used to be (or update these instead of deleting)
    // - add new edges between the new vertices

#ifdef PROFILING
    time->start_section(idx_time_cycledec);
#endif

    std::vector<pdedge> cross_e_;
    for (int cnt = 0, i = 0; cnt < (int)cross_e.size(); ++cnt) {
      cross_e_.push_back(cross_e[i]);
#ifdef DEBUG_CHECKS
      if (!start_e.count(cross_e[i].pj))
        throw std::runtime_error("cycle does not close: next element missing");
#endif
      i = start_e[cross_e[i].pj];
#ifdef DEBUG_CHECKS
      if (cnt == (int)cross_e.size() - 1) {
        if (i != 0)
          throw std::runtime_error("cycle does not close: last element does not lead back to beginning");
      } else {
        if (i == 0)
          throw std::runtime_error("cycle closes too fast");
      }
#endif
    }
    cross_e = cross_e_;

#ifdef PROFILING
    time->start_section(idx_time_add);
#endif

    // add new edges
    pdedge first_e, prev_e;
    bool first = true;
    for (auto e : cross_e) {
      // calculate intersect of e and new halfspace
      // i.e. the segment [mesh.dvert[e.di], mesh.dvert[e.dj]]
      // with the plane with normal vector mesh.pvert[pi] (== hs)
      Vec4 dv;

      auto &a = mesh.pvert[e.pi];
      auto &b = mesh.pvert[e.pj];
      auto &c = mesh.pvert[pi];

      T dt = a({0,1,2}).cross(b({0,1,2})).dot(c({0,1,2}));

      if (std::abs(dt) < 1e-5) {
        // if the determinant is small, the calculation of the intersect from 3 halfplanes would be
        // imprecise, so do an edge-plane intersection instead (this is not always used as it has
        // it's own imprecision problems for very long edges like the ones on the boundary)

        auto &a2 = mesh.dvert[e.di];
        auto &b2 = mesh.dvert[e.dj];

        // <a + l * (b - a), hs> = 0
        T l = std::clamp(hs.dot(a2) / hs.dot(a2 - b2), (T)0., (T)1.);

        dv = a2 + l * (b2 - a2);
      } else {
        dv = Vec4(-a({1,2,3}).cross(b({1,2,3})).dot(c({1,2,3})) / dt,
                  a({2,3,0}).cross(b({2,3,0})).dot(c({2,3,0})) / dt,
                  -a({3,0,1}).cross(b({3,0,1})).dot(c({3,0,1})) / dt,
                  1.0);
      }

      {
        // run some iterations of cyclic projections to refine the result
        int max_its = 100;
        T max_err = 1e-9;

        auto na = a.normalized();
        auto nb = b.normalized();
        auto nc = c.normalized();

        T err = std::numeric_limits<T>::infinity(), d;
        for (int i = 0; i < max_its && err > max_err; ++i) {
          err = 0;
          d = na.dot(dv); dv -= na * d;  err = std::max(err, d);
          d = nb.dot(dv); dv -= nb * d;  err = std::max(err, d);
          d = nc.dot(dv); dv -= nc * d;  err = std::max(err, d);
          dv /= dv(3);
        }
      }

      // replace dual start of crossing edge with new vertex
      e.di = mesh.add_dual(dv);

      // add edge from new vertex to remaining old part
      mesh.add_edge(e);

      if (first) {
        // remember first edge to close the cycle later
        first_e = e;
        first = false;
      } else {
        // if not first edge, connect to previous
        mesh.add_edge({e.pi, pi, prev_e.di, e.di});
      }

      // remember last edge
      prev_e = e;
    }

    if (!first) {
      // close cycle
      mesh.add_edge({first_e.pi, pi, prev_e.di, first_e.di});
    }

#ifdef PROFILING
    time->start_section(idx_time_remove);
#endif

    // remove old cut off dual vertices (including all connecting edges)
    for (const int &di : remove_di)
      mesh.remove_dual(di);

#ifdef PROFILING
    time->end_section();
    remove_count += remove_di.size();
#endif

    return phint;
  }

  /*! Add python bindings to module m.  */
  static void bind(py::module_ &m) {
    py::class_<halfspace_intersection<T>>(m, ("HalfspaceIntersection_" + type_name<T>::value()).c_str())
      .def(py::init<const Eigen::Ref<const halfspace_intersection<T>::Vec3>&,
           const Eigen::Ref<const halfspace_intersection<T>::Vec3>&>(),
           py::arg("lb"), py::arg("ub"))
      .def(py::init<const Eigen::Ref<const halfspace_intersection<T>::Vec3>&,
           const Eigen::Ref<const halfspace_intersection<T>::Vec3>&,
           const Eigen::Ref<const Eigen::Matrix<T, Eigen::Dynamic, 4, Eigen::RowMajor>>>(),
           py::arg("lb"), py::arg("ub"), py::arg("hss"))
      .def("add_halfspace", &halfspace_intersection<T>::add_halfspace,
           py::arg("hs"), py::arg("hint")=-1)
      .def("add_halfspaces", &halfspace_intersection<T>::add_halfspaces,
           py::arg("hss"))
#ifdef PROFILING
      .def_readonly("time", &halfspace_intersection<T>::time)
      .def_readonly("remove_count", &halfspace_intersection<T>::remove_count)
#endif
      .def_readonly("mesh", &halfspace_intersection<T>::mesh);
  }
};

#endif
