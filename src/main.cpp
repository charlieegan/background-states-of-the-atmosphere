#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <list>
#include <algorithm>
#include <numeric>

#if __has_include(<format>)
#include <format>
#define FORMAT std::format
#else
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#define FORMAT fmt::format
#endif

#define DEBUG_CHECKS

#include "laguerre_diagram.hpp"
#include "halfspace_intersection.hpp"


namespace py = pybind11;

class laguerre_diagram_old
{
public:
  double sqr(const double &v) {
    return v * v;
  }
  
protected:
  int standardize(const int &i) {
    if (i < cnt)
      return i;
    else if (i < cnt + 6)
      return cnt;
    else
      return cnt + 1; 
  }

  
public:
  /* 
   * primal indices index halfplanes, index 0 to cnt-1 are for seeds, cnt to cnt+5 are boundaries, cnt+6 to the end are 0-cells
   * dual indices index intersections (corresponding to dual facets)
   * 
   */
  
  typedef Eigen::Ref<const Eigen::Matrix<int,    Eigen::Dynamic, 3, Eigen::RowMajor>> facets_t;
  typedef Eigen::Ref<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> points_t;
  typedef Eigen::Ref<const Eigen::Matrix<int,    Eigen::Dynamic, 2, Eigen::RowMajor>> seeds_t;

  facets_t facets;                          // dual facets (by dual indices)
  points_t intersections;                   // intersections (by dual indices)
  seeds_t seeds;                            // seeds (by primal indices)
  int cnt, dcnt;                            // primal count and dual count
  physical_parameters phys;                 // physical parameters

  std::vector<pdedge> edges;                // primal-dual edges (pair of primal and dual indices belonging together)
  std::vector<std::vector<int>> pfaces;     // primal faces (for each primal index, list edge indices)
  
  laguerre_diagram_old(facets_t facets,         // (dcnt,3) index matrix of dual facets
                       points_t intersections,  // (dcnt,3) position matrix of corresponding intersections
                       seeds_t seeds,           // (cnt,2) seed positions (y1, y2); N < n
                       physical_parameters phys) :
    facets(facets), intersections(intersections), seeds(seeds),
    cnt(seeds.rows()), dcnt(facets.rows()),
    phys(phys), pfaces(cnt) {

    // check input dimensions
    if (dcnt != intersections.rows()) {
      throw std::runtime_error("facets and intersections do not have the same length");
    }

    get_faces_and_edges();
    // calculate_curved_faces();
  }

  std::string repr() const {
    std::stringstream ss;
    ss << "LaguerreDiagram(face_count = " << cnt << ", vertex_count = " << dcnt << ")";
    return ss.str();
  }
  
  void get_faces_and_edges() {
    // for a sorted pair of primal indices, map to dual indices which are adjacent to both of them (to build edges)
    std::map<std::pair<int, int>, int> adj;

    // get primal faces and edges from facets
    for (int di = 0; di < dcnt; di++) {
      for (int j = 0; j < 3; j++) {
        int j2 = (j == 2 ? 0 : j + 1);

        int pi = facets(di, j);
        int pi2 = facets(di, j2);
          
        auto it = adj.find(std::minmax({pi, pi2}));
        if (it != adj.end()) {
          // there is an edge between primal vertices pi and pi2 and dual vertices di and it->second
          
          // if at least one of the primal vertices corresponds to a seed position...
          if (pi < cnt || pi2 < cnt) {
            // store edge
            edges.push_back(pdedge({standardize(pi), standardize(pi2), di, it->second}));

            // and add it to faces
            if (pi < cnt)
              pfaces[pi].push_back(edges.size() - 1);
            if (pi2 < cnt)
              pfaces[pi2].push_back(edges.size() - 1);
          }

          // erase map entry as it is not needed anymore (smaller map should be more efficient)
          adj.erase(it);
        }  else {
          adj[std::minmax({pi, pi2})] = di;
        }
      }
    }

    // fix edge directions (required if input convex hull normals are not unified)
    int pnonempty = -1;
    for (int pi = 0; pi < cnt; ++pi) {
      if (pfaces[pi].size()) {
        pnonempty = pi;
        break;
      }
    }
    
    if (pnonempty == -1)
      throw std::runtime_error("all faces are empty");

    // orient pfaces and edges by traversing them via dfs
    std::vector<int> st({pnonempty});
    std::vector<bool> pused(cnt);
    std::vector<bool> eused(edges.size()); // whether edge was visited (if true it is guaranteed to be correctly oriented)

#ifdef DEBUG_CHECKS
    int itc = 0;
#endif
    
    while (!st.empty()) {
      // take face from stack, if visited continue, otherwise mark visited
      int pi = st.back();

#ifdef DEBUG_CHECKS
      if (pi < 0 || pi >= cnt)
        throw std::runtime_error(FORMAT("dfs tried to visit invalid face {}", pi));
#endif
      
      st.pop_back();
      if (pused[pi])
        continue;
      pused[pi] = true;

       // index of edge to take as orientation reference (some previously visited edge if exists)
      int istart = -1;

      // make bidirectional adjacency map of edges of face and try to find istart
      std::map<int, std::vector<int>> fw;
      for (const int &ei : pfaces[pi]) {
        auto &e = edges[ei];

        fw[e.di].push_back(ei);
        fw[e.dj].push_back(ei);
        if (istart == -1 && eused[ei])
          istart = ei;
      }


      // if no edge was visited previously, just take first edge as reference (should only happen on first face)
      if (istart == -1) {
        istart = pfaces[pi][0];
        eused[istart] = true;
      }

#ifdef DEBUG_CHECKS
      for (auto &[p, eiv] : fw)
        if (eiv.size() != 2)
          throw std::runtime_error(FORMAT("face edges are not closed for face {}", pi));
#endif

      // traverse circle of edges and orient them as necessary
      std::vector<int> nf;
      nf.reserve(pfaces[pi].size());
      for (int ei = istart, i = 0; i < (int)pfaces[pi].size(); ++i) {
        nf.push_back(ei);
        
        // find next edge to visit
        auto &nx = fw[edges[ei].dend(pi)];
        int ej = (nx[0] == ei ? nx[1] : nx[0]);

        // fix orientation of ej if necessary
        if (edges[ej].dstart(pi) != edges[ei].dend(pi)) {
#ifdef DEBUG_CHECKS
          if (eused[ej])
            throw std::runtime_error(FORMAT("edge orientation contradiction at face {} with edges {} and {}", pi, ei, ej));
#endif    
          std::swap(edges[ej].di, edges[ej].dj);
        }

        // mark ej visited
        eused[ej] = true;

        // go to ej
        ei = ej;

#ifdef DEBUG_CHECKS
        if (i == (int)pfaces[pi].size() - 1 && ej != istart)
          throw std::runtime_error(FORMAT("edge cycle did not close at face {}", pi));
        if (i < (int)pfaces[pi].size() - 1 && ej == istart)
          throw std::runtime_error(FORMAT("edge cycle closed too early at face {}", pi));
#endif
      }
      // replace face with ordered version
      pfaces[pi] = nf;

      // go over neighboring faces and add them to stack if not visited already
      for (auto &ei : pfaces[pi]) {
        int pj = (edges[ei].pi == pi ? edges[ei].pj : edges[ei].pi);
        if (pj < cnt && !pused[pj])
          st.push_back(pj);
      }

      
#ifdef DEBUG_CHECKS
      if (itc++ >= cnt)
        throw std::runtime_error("dfs does not terminate correctly");
#endif
    }
    
  }

  void calculate_curved_faces() {
    
  }
};


void hello() {
  py::print("Hello from C++!"); 
}
 

PYBIND11_MODULE(_atmosphere_bgs, m) {
  m.def("_hello", &hello); 

  BIND_PHYSICAL_PARAMETERS(m);
  BIND_SIMULATION_PARAMETERS(m);
  BIND_TANGENT_POINT(m);
  BIND_DISCRETIZED_LINE_SEGMENT(m); 

  BIND_PDEDGE(m);
  BIND_PDMESH(m);
  BIND_HALFSPACE_INTERSECTION(m);

  BIND_DIAGRAM_EDGE(m);
  BIND_LAGUERRE_DIAGRAM(m);
  
  py::class_<laguerre_diagram_old> (m, "_LaguerreDiagram")
    .def(py::init<laguerre_diagram_old::facets_t,
                  laguerre_diagram_old::points_t,
                  laguerre_diagram_old::seeds_t,
                  physical_parameters>(),
         py::arg("dual_facets"),
         py::arg("intersections"),
         py::arg("seeds"),
         py::arg("physical_parameters"))
    .def("__repr__", &laguerre_diagram_old::repr)
    .def_readonly("facets", &laguerre_diagram_old::facets)
    .def_readonly("intersections", &laguerre_diagram_old::intersections)
    .def_readonly("seeds", &laguerre_diagram_old::seeds)
    .def_readonly("cnt", &laguerre_diagram_old::cnt)
    .def_readonly("dcnt", &laguerre_diagram_old::dcnt)
    .def_readonly("physical_parameters", &laguerre_diagram_old::phys)
    .def_readonly("pfaces", &laguerre_diagram_old::pfaces)
    .def_readonly("edges", &laguerre_diagram_old::edges);

    
}
