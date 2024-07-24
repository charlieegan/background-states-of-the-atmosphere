#ifndef HALFSPACE_INTERSECTION_HPP
#define HALFSPACE_INTERSECTION_HPP

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <deque>
#include <Eigen/Dense>
#include "dvector.hpp"

// an (oriented) primal-dual edge represented by primal and dual indices
struct pdedge {
  int pi, pj; // primal indices neighboring
  int di, dj; // corresponding dual indices neighboring

  // pdedge(const pdedge&) = default;
  // pdedge(pdedge&&) = default;
  // pdedge& operator=(const pdedge&) = default;
  // pdedge& operator=(pdedge&&) = default;
  
  // get start of dual edge relative to primal index p in [pi, pj]
  const int &dstart(const int &p) const {
#ifdef DEBUG_CHECKS
    if (p != pi && p != pj)
      throw std::runtime_error("dstart called with primal index that is not part of edge");
#endif
    return (p == pi) ? di : dj;
  }
  // get end of dual edge relative to primal index p in [pi, pj]
  const int &dend(const int &p) const {
#ifdef DEBUG_CHECKS
    if (p != pi && p != pj)
      throw std::runtime_error("dend called with primal index that is not part of edge");
#endif
    return (p == pi) ? dj : di;
  }
  // get start of primal edge relative to dual index d in [di, dj]
  const int &pstart(const int &d) const {
#ifdef DEBUG_CHECKS
    if (d != di && d != dj)
      throw std::runtime_error("pstart called with dual index that is not part of edge");
#endif
    return (d == di) ? pi : pj;
  }
  // get end of primal edge relative to dual index d in [di, dj]
  const int &pend(const int &d) const {
#ifdef DEBUG_CHECKS
    if (d != di && d != dj)
      throw std::runtime_error("pend called with dual index that is not part of edge");
#endif
    return (d == di) ? pj : pi;
  }
  
  std::string repr() const {
    std::stringstream ss;
    ss << "PDEdge(" << pi << ", " << pj << "; " << di << ", " << dj << ")";
    return ss.str();
  }

  pdedge mirrored() const {
    return {pj, pi, dj, di};
  }
  
  struct cmp_pj { bool operator()(const pdedge &e, const pdedge &f) const { return e.pj < f.pj; } };
  struct cmp_dj { bool operator()(const pdedge &e, const pdedge &f) const { return e.dj < f.dj; } };


  // order the given edges to make up a cycle on primal side (if primal == true) or dual side (if primal == false)
  // if the edges do not make up a cycle, a runtime_error will the thrown
  template <bool primal = true>
  static std::vector<pdedge> order_cycle(const std::vector<pdedge> &edg) {
    if (edg.size() <= 2)
      return edg;
    std::map<int, int> next;
    for (int i = 0; i < (int)edg.size(); ++i)
      next[primal ? edg[i].pi : edg[i].di] = i;
    std::vector<pdedge> res;
    res.reserve(edg.size());
    res.push_back(edg[0]);
    for (int i = 1, j = 0; i < (int)edg.size(); ++i) {
      auto it = next.find(primal ? edg[j].pj : edg[j].dj);
      if (it == next.end()) {
        std::stringstream ss;
        ss << "edges do not make up a cycle: ";
        for (auto & e : edg)
          ss << e.repr() << " ";
        throw std::runtime_error(ss.str());
      }
      j = it->second;
      res.push_back(edg[j]);
    }
    if ((primal ? res.back().pj : res.back().dj) != (primal ? res.front().pi : res.front().di)) {
      throw std::runtime_error("edges do not make up a cycle (cycle does not close)");
    }
    return res;
  }
};

#define BIND_PDEDGE(m)                                  \
  py::class_<pdedge>(m, "PDEdge")                       \
  .def("__repr__", &pdedge::repr)                       \
  .def("pstart", &pdedge::pstart)                       \
  .def("pend", &pdedge::pend)                           \
  .def("dstart", &pdedge::dstart)                       \
  .def("dend", &pdedge::dend)                           \
  .def_readonly("pi", &pdedge::pi)                      \
  .def_readonly("pj", &pdedge::pj)                      \
  .def_readonly("di", &pdedge::di)                      \
  .def_readonly("dj", &pdedge::dj);


struct pdmesh {
  // primal adjacency list - for each primal index: list of incident edges
  std::vector<std::vector<pdedge>> padj;

  // dual adjacency list - for each dual index: list of incident edges (3 for a dual-triangle-mesh)
  dvector<std::vector<pdedge>> dadj;

  // primal vertices
  std::vector<Eigen::Vector4d> pvert;
  
  // dual vertices
  dvector<Eigen::Vector4d> dvert;
  
  pdmesh(const int &pcnt = 0, const int &dcnt = 0) :
    padj(pcnt), dadj(dcnt),
    pvert(pcnt), dvert(dcnt) {}

  static pdmesh cube(const Eigen::Ref<const Eigen::Vector3d> &lb,
                     const Eigen::Ref<const Eigen::Vector3d> &ub) {
    pdmesh res(6, 8);

    // a cuboid mesh has 6 prima vertices, 8 dual vertices, 12 edges

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

    return res;
  }
  
  int pcnt() const {
#ifdef DEBUG_CHECKS
    if (padj.size() != pvert.size())
      throw std::runtime_error("padj and pvert do not have the same size");
#endif
    return padj.size();
  }

  int dcnt() const {
#ifdef DEBUG_CHECKS
    if (dadj.size() != dvert.size())
      throw std::runtime_error("dadj and dvert do not have the same size");
#endif
    return dadj.size();
  }
  
  // add primal-dual edge to mesh
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

  // remove primal-dual edge from mesh 
  void remove_edge(const pdedge &e) {
    // TODO: efficiency (using vectors may actually be best when number of elements is very small -> benchmark)
#ifdef DEBUG_CHECKS
    int i0 = 
#endif
    std::erase_if(padj[e.pi], [&e](const pdedge &f){ return f.pj == e.pj; });
#ifdef DEBUG_CHECKS
    int i1 = 
#endif
    std::erase_if(padj[e.pj], [&e](const pdedge &f){ return f.pj == e.pi; });
#ifdef DEBUG_CHECKS
    int i2 = 
#endif
    std::erase_if(dadj[e.di], [&e](const pdedge &f){ return f.dj == e.dj; });
#ifdef DEBUG_CHECKS
    int i3 = 
#endif
    std::erase_if(dadj[e.dj], [&e](const pdedge &f){ return f.dj == e.di; });

#ifdef DEBUG_CHECKS
    if (i0 < 1 || i1 < 1 || i2 < 1 || i3 < 1)
      throw std::runtime_error("tried to erase edge that does not exist");
    if (i0 > 1 || i1 > 1 || i2 > 1 || i3 > 1)
      throw std::runtime_error("erase removed multiple edges");
#endif
  }

  // remove dual vertex at index with all corresponding edges
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

  // add a new dual, return the index
  int add_dual() {
    dvert.add(Eigen::Vector4d());
    int di = dadj.add(std::vector<pdedge>());
    // std::cout << "added dual " << di << std::endl;
    return di;
  }

  // ad a new primal vertex and return the index
  int add_primal(const Eigen::Ref<const Eigen::Vector4d> &hs) {
    pvert.push_back(hs);
    padj.push_back(std::vector<pdedge>());
    return (int)pvert.size() - 1;
  }
  
  // make sure the order of dadj[di] is correct if dadj[di].size() == 3
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
      throw std::runtime_error("dual edges to not form cycle");
#endif

    auto &a = pvert[dadj[di][0].pi];
    auto &b = pvert[dadj[di][1].pi];
    auto &c = pvert[dadj[di][2].pi];         
    
    dvert[di] =
      {
        -a({1,2,3}).cross(b({1,2,3})).dot(c({1,2,3})),
        a({2,3,0}).cross(b({2,3,0})).dot(c({2,3,0})),
        -a({3,0,1}).cross(b({3,0,1})).dot(c({3,0,1})),
        a({0,1,2}).cross(b({0,1,2})).dot(c({0,1,2}))
      };

  }


  // find how far inside / outside of the halfspace corresponding to pi the dual point di is
  // negative if inside, positive outside
  double hs_error(const int &pi, const int &di) const {
    return pvert[pi].dot(dvert[di]) / (dvert[di](3) * pvert[pi].head<3>().norm());
  }
  
  // check whether halfspace defined by pi contains dual point di
  // if di is not a valid dual index, return default_return
  bool contains(const int &pi, const int &di) const {
#ifdef DEBUG_CHECKS
    if (pi < 0 || pi >= pcnt() || !dadj.contains_idx(di))
      throw std::runtime_error("called contains() with invalid points");
#endif
    return pvert[pi].dot(dvert[di]) <= 0;
  }
  
  // get neighbors of primal index pi (as edges e with e.pi == pi)
  const std::vector<pdedge>& pneigh(const int &pi) const {
#ifdef DEBUG_CHECKS
    if (pi < 0 || pi >= pcnt())
      throw std::runtime_error("called pneigh with invalid index");
#endif
    return padj[pi];
  }

  // get neighbors of dual index di (as edges e with e.di == di)
  const std::vector<pdedge>& dneigh(const int &di) const {
#ifdef DEBUG_CHECKS
    if (!dadj.contains_idx(di))
      throw std::runtime_error("called pneigh with invalid index");
#endif
    return dadj[di];
  }

  std::string repr() const {
    std::stringstream ss;
    ss << "PDMesh(pcnt = " << pcnt() << ", dcnt = " << dcnt() << ")";
    return ss.str();
  }
  
  std::vector<Eigen::Vector4d> primal_poly(const int &pi) const {
    std::vector<Eigen::Vector4d> res;
    for (auto &e : pdedge::order_cycle<false>(padj[pi])) {
      res.push_back(dvert[e.di]);
    }
    return res;
  }
  
  bool check_integrity() const {
    bool res = true;

    if (dadj.size() != dvert.size()) {
      std::cout << "dadj and dvert have different sizes" << std::endl;
      res = false;
    }
    if (padj.size() != pvert.size()) {
      std::cout << "padj and pvert have different sizes" << std::endl;
      res = false;                                                     
    }
    
    for (const auto &v : dadj) {
      // check dual face size
      if (v.size() != 3) {
        std::cout << "non-triangular dual face found" << std::endl;
        res = false;
      }

      // check dual face orientation / order
      int first_pi, prev_pj;
      bool first = true;
      for (const auto &e : v) {
        if (first) {
          first_pi = e.pi;
          first = false;
        } else {
          if (e.pi != prev_pj) {
            std::cout << "incorrectly ordered or oriented dual face found" << std::endl;
            res = false;
          }
        }
        prev_pj = e.pj;
      }
      if (first_pi != prev_pj) {
        std::cout << "incorrectly ordered dual face found" << std::endl;
        res = false;
      }
    }
    for (const auto &v : padj) {
      // check primal face size
      if (v.size() == 1 || v.size() == 2) {
        std::cout << "degenerated primal face found" << std::endl;
        res = false;
      }
    }
   
    return res;
  }
};

#define BIND_PDMESH(m)                                                  \
  py::class_<pdmesh>(m, "PDMesh")                                       \
  .def("__repr__", &pdmesh::repr)                                       \
  .def("check_integrity", &pdmesh::check_integrity)                     \
  .def("primal_poly", &pdmesh::primal_poly)                             \
  .def_property_readonly("pcnt", &pdmesh::pcnt)                         \
  .def_property_readonly("dcnt", &pdmesh::dcnt)                         \
  .def_readonly("padj", &pdmesh::padj)                                  \
  .def_property_readonly("dadj",                                        \
                         [](const pdmesh *m){ return m->dadj.to_vector();}) \
  .def_readonly("pvert", &pdmesh::pvert)                                \
  .def_property_readonly("dvert",                                       \
                         [](const pdmesh *m){ return m->dvert.to_vector();});


class halfspace_intersection
{
public:

  pdmesh mesh;
  
  std::vector<int> used;
  int used_ctr;
  

  // start with cuboid domain defined by lower and upper bounds
  halfspace_intersection(const Eigen::Ref<const Eigen::Vector3d> &lb = Eigen::Vector3d(0,0,0),
                         const Eigen::Ref<const Eigen::Vector3d> &ub = Eigen::Vector3d(1,1,1)) :
    mesh(pdmesh::cube(lb, ub)), used(8), used_ctr(0) { }

  halfspace_intersection(const Eigen::Ref<const Eigen::Vector3d> &lb,
                         const Eigen::Ref<const Eigen::Vector3d> &ub,
                         const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor>> &hss) :
    halfspace_intersection(lb, ub) {
    add_halfspaces(hss);
  }
  
  void add_halfspaces(const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor>> &hss) {
    for (auto &hs : hss.rowwise())
      add_halfspace(hs);
  }
  
  void add_halfspace(const Eigen::Ref<const Eigen::Vector4d> &hs) {
    // increment used_ctr to effectively get fresh used array
    used_ctr++;
    // reset used array if necessary (should not really happen with reasonable input sizes)
    if (used_ctr < 0) {
      std::fill(used.begin(), used.end(), 0);
      used_ctr = 1;
    }
    // make sure used array is large enough
    while ((int)used.size() < mesh.dadj.capacity()) used.push_back(0);

    // add new primal
    int pi = mesh.add_primal(hs);

    // dual vertex that is not contained in hs
    int di_start = -1;

    // find dual point outside of halfspace
#if true
    int di = *mesh.dadj.begin_idx(), steps=0;
    double err = mesh.hs_error(pi, di);
    while (true) {
      ++steps;

      used[di] = used_ctr;
      
      if (err >= 0) {
        di_start = di;
        break;
      }
      
      double max_err = -std::numeric_limits<double>::infinity(), nerr;
      int ndi = -1;
      for (const auto &e : mesh.dneigh(di))
        if ((nerr = mesh.hs_error(pi, e.dj)) > max_err && used[e.dj] != used_ctr)
          max_err = nerr, ndi = e.dj;

      // TODO: this can break due to precision problems -> fix
      if (std::isinf(max_err)) {
        // py::print("break: new error", max_err, "cur error", err);
        di_start = -1;
        break;
      }
      
      di = ndi;
      err = max_err;
    }
    // py::print("for primal", pi,
    //           "found dual index", di_start,
    //           "with err", err,
    //           "contains:", (di_start >= 0 ? mesh.contains(pi, di_start) : false),
    //           "in", steps, "steps");
    // std::cout << "found dual index " << di << " with err " << err << " contains: " << (di >= 0 ? mesh.contains(pi, di) : false) << std::endl;

#else
    for (auto it = mesh.dadj.begin_idx(); it != mesh.dadj.end_idx(); ++it) {
      int di = *it;
      // std::cout << "check contains(" << pi << ", " << di << ")" << std::endl;
      if (!mesh.contains(pi, di)) {
        di_start = di;
        break;
      }
    }
#endif

    // if there is no point outside, the new halfspace does nothing
    if (di_start < 0) {
      // std::cout << "halfspace completely contains domain" << std::endl;
      return;
    }

    used_ctr++;
    used[di_start] = used_ctr;
    std::vector<int> remove_di({di_start});  // dual vertices that are cut off
    std::vector<pdedge> cross_e;             // edges crossing from cut off part to remaining part (e.di is removed, e.dj is remaining)

    // (bfs) traverse set of cut off vertices and store the vertices and edges connecting them to the remaining part
    std::deque<int> q({di_start});
    while (!q.empty()) {
      int di = q.front();
      q.pop_front();

      // std::cout << "visit " << di << std::endl;

      for (const auto &e : mesh.dneigh(di)) {
        // ignore already visited vertices
        if (used[e.dj] == used_ctr)
          continue;
        
        if (!mesh.contains(pi, e.dj)) {
          // e.dj is also cut off -> store and push onto queue
          q.push_back(e.dj);
          remove_di.push_back(e.dj);
          // mark used
          used[e.dj] = used_ctr;
        } else {
          // e.dj is not cut off -> store the edge
          cross_e.push_back(e);
        }

      }
    }
    
    // std::cout << "cut verts: ";
    // for (auto &x: remove_di)
    //   std::cout << x << " ";
    // std::cout << std::endl << "cut edges: ";
    // for (auto &e : cross_e)
    //   std::cout << e.repr() << " ";
    // std::cout << std::endl;

    // - remove cut dual vertices with all their edges
    // - add new dual vertices at intersections
    // - add new edges where cross_e used to be (or update these instead of deleting)
    // - add new edges between the new vertices

    // remove cut off dual vertices (including all connecting edges)
    for (const int &di : remove_di)
      mesh.remove_dual(di);

    // order the cut edges by primal indices (i.e. such that e[i].pj == e[i+1].pi)
    cross_e = pdedge::order_cycle<true>(cross_e);

    // add new edges
    pdedge first_e, prev_e;
    bool first = true;
    for (auto &e : cross_e) {
      // replace dual start of crossing edge with new vertex
      e.di = mesh.add_dual();

      // add edge from new vertex to remaining old part
      // std::cout << "add edge " << e.repr() << " between new and old" << std::endl;
      mesh.add_edge(e);
      
      if (first) {
        // remember first edge to close the cycle later
        first_e = e;
        first = false;
      } else {
        // if not first edge, connect to previous
        // std::cout << "add edge between new and new" << std::endl;
        mesh.add_edge({e.pi, pi, prev_e.di, e.di});
      }

      // remember last edge
      prev_e = e;
    }

    // close cycle
    if (!first) {
      // std::cout << "add closing edge" << std::endl;
      mesh.add_edge({first_e.pi, pi, prev_e.di, first_e.di});
    }
  }
};

#define BIND_HALFSPACE_INTERSECTION(m)                                  \
  py::class_<halfspace_intersection>(m, "HalfspaceIntersection")        \
  .def(py::init<const Eigen::Ref<const Eigen::Vector3d>&,               \
                const Eigen::Ref<const Eigen::Vector3d>&>(),            \
       py::arg("lb"), py::arg("ub"))                                    \
  .def(py::init<const Eigen::Ref<const Eigen::Vector3d>&,               \
                const Eigen::Ref<const Eigen::Vector3d>&,               \
                const Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 4, Eigen::RowMajor>>>(), \
       py::arg("lb"), py::arg("ub"), py::arg("hss"))                    \
  .def("add_halfspace", &halfspace_intersection::add_halfspace,         \
       py::arg("hs"))                                                   \
  .def("add_halfspaces", &halfspace_intersection::add_halfspaces,       \
       py::arg("hss"))                                                  \
  .def_readonly("mesh", &halfspace_intersection::mesh);




#endif

