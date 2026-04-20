#define PROFILING

#include "laguerre_diagram.hpp"

// explicit instantiations
template class laguerre_diagram<double>;
template class laguerre_diagram<long double>;

// ################################# diagram_edge ##################################
template <typename T>
std::string laguerre_diagram<T>::diagram_edge::repr() const {
  std::stringstream ss;
  ss << "diagram_edge("
     << "pi = " << pi
     << ", pj = " << pj
     << ", di = " << di
     << ", dj = " << dj
     << ")";
  return ss.str();
}

template <typename T>
void laguerre_diagram<T>::diagram_edge::bind(py::module_ &m) {
  py::class_<laguerre_diagram<T>::diagram_edge>(m, ("DiagramEdge_" + type_name<T>::value()).c_str())
    .def("__repr__", &laguerre_diagram<T>::diagram_edge::repr)
    .def_readonly("pi", &laguerre_diagram<T>::diagram_edge::pi)
    .def_readonly("pj", &laguerre_diagram<T>::diagram_edge::pj)
    .def_readonly("di", &laguerre_diagram<T>::diagram_edge::di)
    .def_readonly("dj", &laguerre_diagram<T>::diagram_edge::dj)
    .def_readonly("ls", &laguerre_diagram<T>::diagram_edge::ls)
    .def_readonly("dif", &laguerre_diagram<T>::diagram_edge::dif);
}


#ifdef PROFILING
template <typename T>
void laguerre_diagram<T>::setup_timer() {
  if (!time)
    time = std::make_shared<timer>();
  idx_time_halfspace_intersection = time->get_index_from_name("laguerre_diagram.halfspace_intersection");
  idx_time_extract_vertices = time->get_index_from_name("laguerre_diagram.extract_vertices");
  idx_time_extract_edges = time->get_index_from_name("laguerre_diagram.extract_edges");
  idx_time_calc_areas = time->get_index_from_name("laguerre_diagram.calc_areas");
  idx_time_calc_integrals = time->get_index_from_name("laguerre_diagram.calc_integrals");
  idx_time_jacobian = time->get_index_from_name("laguerre_diagram.jacobian");
  idx_time_prepare_rasterizer = time->get_index_from_name("laguerre_diagram.prepare_rasterizer");
  idx_time_touching_dual = time->get_index_from_name("laguerre_diagram.touching_dual");
  idx_time_touching_dual_n = time->get_index_from_name("laguerre_diagram.touching_dual_n");
}
#endif

template <typename T>
laguerre_diagram<T>::laguerre_diagram(const seeds_t &ys,
                                      const Eigen::Ref<const VecX> &duals,
                                      std::shared_ptr<physical_parameters> phys,
                                      const simulation_parameters &sim) :
  parent(NULL), n(ys.rows()), ys(ys), duals(duals), phys(phys), sim(sim), hints(n + std::max(0, sim.boundary_res), -1) {

#ifdef PROFILING
  setup_timer();
#endif

  if (duals.array().isNaN().any())
    throw std::runtime_error("duals contain NaN");
  
  if ((int)duals.size() != n + 1)
    throw std::runtime_error("the size of the dual vector must be one larger than the first dimension of seed positions (the last entry of duals corresponds to the boundary cell)");

  do_hs_intersect();
  
  extract_diagram();
}

template <typename T>
laguerre_diagram<T>::laguerre_diagram(std::shared_ptr<laguerre_diagram<T>> parent,
                                      const Eigen::Ref<const VecX> &duals,
                                      const bool &keep_parent) :
  parent(keep_parent ? parent : NULL), n(parent->n), ys(parent->ys), duals(duals), phys(parent->phys), sim(parent->sim), hints(parent->hints) {

#ifdef PROFILING
  setup_timer();
#endif

  if (duals.array().isNaN().any())
    throw std::runtime_error("duals contain NaN");

  if ((int)duals.size() != n + 1)
    throw std::runtime_error("the size of the dual vector must be one larger than the first dimension of seed positions (the last entry of duals corresponds to the boundary cell)");

  do_hs_intersect();

  extract_diagram();
}

template <typename T>
void laguerre_diagram<T>::detach() {
  parent = NULL;
}

template <typename T>
void laguerre_diagram<T>::do_hs_intersect() {
#ifdef PROFILING
  time->start_section(idx_time_halfspace_intersection);
#endif

  // perform the halfspace intersection

  Vec3 zeta_min, zeta_max;
  zeta_min << phys->tf<T>(sim.spmin.cast<T>()), (T)-1e100;
  zeta_max << phys->tf<T>(sim.spmax.cast<T>()), (T)1e100;
  hs = halfspace_intersection<T>(zeta_min, zeta_max);

  // add seed halfspaces
  for (int i = 0; i < n; i++) {
    Vec4 H(-0.5 * sqr(ys(i, 0) / phys->a),
           -phys->cp * ys(i, 1),
           1.0,
           phys->Omega * ys(i, 0) + duals(i));

#ifdef DEBUG_CHECKS
    if (i < 0 || i >= hints.size())
      throw std::runtime_error(FORMAT("index in hints oob : {} >= {}", i, hints.size()));
#endif
    
    hints[i] = hs.add_halfspace(H, hints[i]);
  }

  // add boundary halfspaces
  if (sim.boundary_res > 0) {
    boundary_spt.resize(sim.boundary_res);
    
    // evenly spaced halfspaces
    for (int i = 0; i < sim.boundary_res; i++) {
      T s = sim.spmin(0) + (sim.spmax(0) - sim.spmin(0)) * i / std::max(1, sim.boundary_res - 1);
      T z0 = 1. / (1 - s * s);
      boundary_spt[i] = z0;
      
      Vec4 H(-0.5 * sqr(phys->Omega * phys->a / z0), 0.0, 1.0, duals(n) + sqr(phys->Omega * phys->a) / z0);
      
#ifdef DEBUG_CHECKS
      if (n + i < 0 || n + i >= hints.size())
        throw std::runtime_error(FORMAT("index in hints oob : {} >= {}", n + i, hints.size()));
#endif
      
      hints[n + i] = hs.add_halfspace(H, hints[n + i]);
    }
  } else {
    // halfspaces on intersections

    // positions at which to add halfspaces (later)
    T z0gmin = 1. / (1 - sqr(sim.spmin(0)));
    T z0gmax = 1. / (1 - sqr(sim.spmax(0)));
    std::vector<T> z0s({z0gmin, z0gmax});

    // go over all edges
    for (int i = 0; i < n; i++) {
      for (const auto &e : hs.mesh.pneigh(i + 6)) {
        // only consider edges with pi < pj (so each edge only appears once)
        if (e.pi > e.pj)
          continue;

        // do not consider edges on top boundary
        if (e.pi >= n + 6 || e.pj >= n + 6)
          continue;
        
        // only consider inside edges
        if (e.pi < 6 || e.pj < 6)
          continue;
        
        // calculate intersection
        auto yi = ys.row(e.pi - 6);
        auto yj = ys.row(e.pj - 6);
        auto dy = yi - yj;
        auto psii = duals[e.pi - 6];
        auto psij = duals[e.pj - 6];
        auto psitop = duals[n];

        //   z0 * (yi[0]**2 / (2 * a**2)) + z1 * cp * yi[1] - Omega * yi[0] - psii
        // = z0 * (yj[0]**2 / (2 * a**2)) + z1 * cp * yj[1] - Omega * yj[0] - psij
        // = -Omega**2 * a**2 / (2 * z0) - psitop

        auto [z1min, z1max] = std::minmax(hs.mesh.dvert[e.di][1], hs.mesh.dvert[e.dj][1]);
        auto [z0min, z0max] = std::minmax(hs.mesh.dvert[e.di][0], hs.mesh.dvert[e.dj][0]);

        if (dy[1] == 0 || hs.mesh.dvert[e.di][0] == hs.mesh.dvert[e.dj][0]) {
          // vertical line segment (horizontal y, normals)
          T z0 = hs.mesh.dvert[e.di][0]; // 2 * sqr(phys->a) * (phys->Omega * dy[0] + psii - psij) / (sqr(yi[0]) - sqr(yj[0])); // = hs.mesh.dvert[e.di][0]
          T z1 = (-sqr(phys->Omega * phys->a) / (2 * z0)
                  - sqr(yi[0]) / (2 * sqr(phys->a)) * z0
                  + phys->Omega * yi[0] + psii - psitop) / (phys->cp * yi[1]);
          if (z1min <= z1 && z1 <= z1max) {
            z0s.push_back(z0);
          }
        } else {
          // non-vertical line segment
          T c0 = (sqr(yj[0]) * yi[1] - sqr(yi[0]) * yj[1]) / (sqr(phys->a) * dy[1]);
          T c1 = (yi[1] * (phys->Omega * dy[0] + psii - psij) + dy[1] * (-phys->Omega * yi[0] + psitop - psii)) / dy[1];
          T c2 = sqr(phys->Omega * phys->a);
          // intersections are at solutions to c0 * z0**2 + 2 * c1 * z0 + c2 = 0

          // solutions are at z0 = p +- sqrt(q)
          T p = -c1 / c0;
          T q = sqr(p) - c2 / c0;

          if (q > 0) {
            // two intersection candidates - check and add
            q = std::sqrt(q);
            T z0 = p - q;
            if (z0min <= z0 && z0 <= z0max) {
              T z1 = (-sqr(phys->Omega * phys->a) / (2 * z0)
                      - sqr(yi[0]) / (2 * sqr(phys->a)) * z0
                      + phys->Omega * yi[0] + psii - psitop) / (phys->cp * yi[1]);
              // T z1 = (phys->Omega * dy[0] + psii - psij - z0 * (sqr(yi) - sqr(yj)) / (2 * phys->a)) / (phys->cp * dy[1]);
              if (z1min <= z1 && z1 <= z1max) {
                z0s.push_back(z0);
              }
            }
            z0 = p + q;
            if (z0min <= z0 && z0 <= z0max) {
              T z1 = (-sqr(phys->Omega * phys->a) / (2 * z0)
                      - sqr(yi[0]) / (2 * sqr(phys->a)) * z0
                      + phys->Omega * yi[0] + psii - psitop) / (phys->cp * yi[1]);
              if (z1min <= z1 && z1 <= z1max) {
                z0s.push_back(z0);
              }
            }
          } else if (q == 0) {
            // one intersection candidates - check and add
            T z0 = p;
            if (z0min <= z0 && z0 <= z0max) {
              T z1 = (-sqr(phys->Omega * phys->a) / (2 * z0)
                      - sqr(yi[0]) / (2 * sqr(phys->a)) * z0
                      + phys->Omega * yi[0] + psii - psitop) / (phys->cp * yi[1]);
              if (z1min <= z1 && z1 <= z1max) {
                z0s.push_back(z0);
              }
            }
          } else {
            // no intersections - nothing to do
          }
        }
      }
    }

    // sort cells and remove duplicates (just in case)
    std::sort(z0s.begin(), z0s.end());
    int unique = std::unique(z0s.begin(), z0s.end()) - z0s.begin();
    // py::print(FORMAT("found {} intersections of which {} are unique", z0s.size(), unique));
    z0s.resize(unique);

    boundary_spt.resize(z0s.size());
    for (int i = 0; i < (int)z0s.size(); ++i)
      boundary_spt[i] = z0s[i];
    
    // add halfspaces (use previous as hint as they are neighboring, first hs is added without hint)
    if (sim.boundary_res != 0) {
      int hint = -1;
      for (auto &z0 : z0s) {
        Vec4 H(-0.5 * sqr(phys->Omega * phys->a / z0), 0.0, 1.0, duals(n) + sqr(phys->Omega * phys->a) / z0);
        hint = hs.add_halfspace(H, hint);
      }
    }

  }

#ifdef PROFILING
  time->end_section();
#endif
}

template <typename T>
void laguerre_diagram<T>::extract_diagram() {
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

        if (sim.boundary_res < 0 && e.pj >= n + 6) {
          // top edge
          double z0 = verts.row(e.di)[0]; // one end first coordinate
          double z1 = verts.row(e.dj)[0]; // other end first coordinate
          auto yi = ys.row(e.pi - 6); // seed position
          double curve = -0.5 * sqr(phys->Omega * phys->a) / (phys->cp * yi[1]); // inverse factor
          double lin = -0.5 * sqr(yi[0] / phys->a) / (phys->cp * yi[1]); // linear factor
          double cons = (phys->Omega * yi[0] + duals[e.pi - 6] - duals[n]) / (phys->cp * yi[1]); // constant offset
          discretized_line_segment<double>::Vector2 start
            = discretized_line_segment<double>::Vector2(z0, curve / z0 + lin * z0 + cons); // start point of segment
          discretized_line_segment<double>::Vector2 end
            = discretized_line_segment<double>::Vector2(z1, curve / z1 + lin * z1 + cons); // end point of segment
          edglist.push_back({e.pi, e.pj, e.di, e.dj,
              discretized_line_segment<double>(start, end, phys, sim, curve), 0});
        } else {
          // inside edge
          edglist.push_back({e.pi, e.pj, e.di, e.dj,
              discretized_line_segment<double>(verts.row(e.di), verts.row(e.dj), phys, sim), 0});
        }

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
    for (const auto &ei : faces[i + 6]) {
      const auto &e = edglist[ei];
      areas(i) += (e.pi == i + 6 ? e.ls.area : -e.ls.area);
      areaerrs(i) += e.ls.errb;
    }

    if (areas(i) == 0 || areaerrs(i) / areas(i) < sim.area_tolerance) {
      // py::print(i, "- error already in tol:", areaerrs(i) / areas(i));
      continue;
    }

    std::multimap<double, int> errmap;

    for (const auto &ei : faces[i + 6]) {
      const auto &e = edglist[ei];
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

  // handle empty / negative areas
  if (sim.negative_area_scaling > 0) {
    for (int i = 0; i < n; ++i) {
      if (areas(i) <= 1e-9) {
        // halfspace normal vector
        auto pp = hs.mesh.pvert[i + 6].head(3);
        // offset = Omega * ys(i, 0) + duals(i)
        auto os = hs.mesh.pvert[i + 6](3);
        // find closest dual point
        int di = hs.mesh.extremal_dual(pp); // TODO: use hint (has to be converted from primal to dual first)
        auto dp = hs.mesh.dvert[di].head(3) / hs.mesh.dvert[di][3];
        // distance (in dual) for halfspace to touch
        auto dist = -pp.dot(dp) - os;

        // py::print(FORMAT("primal {} at dist {} to dual {}", i, dist, di));
        
        if (dist <= 0)
          continue;
        
        // go over primals that are neighboring to closest dual point and add areas
        auto narea = sim.negative_area_scaling * sqr(dist);
        auto dif = 2 * sim.negative_area_scaling * dist;

        // to just add negative area, take mass from top cell, use the 2 lines below:
        // areas(i) -= narea;
        // vedglist.push_back({i + 6, n + 6, dif); 

        // add negative area taken from neighboring cells
        auto w = hs.mesh.dneigh(di).size();
        for (const auto &e : hs.mesh.dneigh(di)) {
          // there is only a change for non-simulation-boundary neighbors
          if (e.pi >= 6) {
            // subtract area from this, add virtual edge
            areas(i) -= w * narea;
            vedglist.push_back({i + 6, e.pi, (double)(w * dif)});

            // if the other cell is a true cell, add the area there
            if (e.pi - 6 < n)
              areas(e.pi - 6) += w * narea;
          }
        }
      }
    }
  }

#ifdef PROFILING
  time->start_section(idx_time_calc_integrals);
#endif

  // calculate line integrals
  for (auto &e : edglist) {

    e.dif = 0;

    auto yi = ys.row(e.pi-6);

    if (e.pj >= 6) {

      Eigen::Vector2d dcp, dc;
      for (int k = 0; k < e.ls.size(); ++k) {
        // d_x (c(x,yi) - c(x,yj))

        if (e.pj - 6 < n) {
          // inner edge
          auto yj = ys.row(e.pj-6);
          dc = Eigen::Vector2d((sqr(yi(0)) - sqr(yj(0))) * e.ls.x(k,0) / (sqr(phys->a) * sqr(1 - sqr(e.ls.x(k,0)))),
                               phys->cp * phys->kappa / phys->p00 * (yi(1) - yj(1)) * std::pow(e.ls.x(k,1) / phys->p00, phys->kappa-1));
        } else {
          // top boundary edge
          dc = Eigen::Vector2d(sqr(yi(0)) * e.ls.x(k,0) / (sqr(phys->a) * sqr(1 - sqr(e.ls.x(k,0)))) - sqr(phys->Omega * phys->a) * e.ls.x(k,0),
                               phys->cp * phys->kappa / phys->p00 * yi(1) * std::pow(e.ls.x(k,1) / phys->p00, phys->kappa-1));
        }
        
        if (k > 0) {
          double xlen = (e.ls.x.row(k) - e.ls.x.row(k - 1)).norm();
          double hlen = (dcp - dc).norm();
          double i1 = dcp.dot(dc - dcp);
          double i2 = dc.dot(dc - dcp);
          double s = (i1 < 0 ? -1.0 : 1.0);

          if (xlen == 0) {
            // if segment has 0 length, don't add anything
          } else if (hlen < 1e-100) {
            // if dc is very close to constant, just approximate it as constant
            e.dif += xlen / dc.norm();
          } else {
            // otherwise solve the line integral
            e.dif += s * (xlen / hlen) * std::log((dc.norm() * hlen + s * i2) / (dcp.norm() * hlen + s * i1));
          }
        }
        
        dcp = dc;
      }
      
    }

#ifdef DEBUG_CHECKS
    if (std::isinf(e.dif) || std::isnan(e.dif))
      throw std::runtime_error(FORMAT("edge {} dif is {}", e.repr(), e.dif));
#endif
    
  }

#ifdef PROFILING
  time->end_section();
#endif
}

template<typename T>
laguerre_diagram<T>::VecX laguerre_diagram<T>::touching_dual(bool randomize) {
  // get the dual that touches the polytope
  // for primal coords that are part of the polytope, this is unchanged
  // for primal coords that are not part, it is larger
  // equivalent to c-transform on other dual

#ifdef PROFILING
  time->start_section(idx_time_touching_dual);
#endif

  VecX res = duals;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.25, 0.75);
  T top_cell_shift = 1e4; // max amount by which to shift top cell when calculating touching value

  for (int pi = 0; pi < n; ++pi) {
    if (hs.mesh.pneigh(pi + 6).size() == 0) {
      // plane does not currently intersect...

      // direction of halfplane at index pi
      auto pp = hs.mesh.pvert[pi + 6].head(3);
      // z-offset of plane (without dual potential part)
      auto os = phys->Omega * ys(pi, 0);
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

#ifdef PROFILING
  time->start_section(idx_time_touching_dual_n);
#endif

  // top boundary dual
  bool touches = false;
  for (int pi = hs.mesh.pcnt()-1; pi >= n + 6; --pi)
    if (hs.mesh.pneigh(pi).size() > 0) {
      touches = true;
      break;
    }

  if (!touches) {
    T mx = std::numeric_limits<T>::infinity(), mn = std::numeric_limits<T>::infinity();
    for (int pi = hs.mesh.pcnt()-1; pi >= n + 6; --pi) {
      auto pp = hs.mesh.pvert[pi].head(3);
      auto os = hs.mesh.pvert[pi][3] - duals[n];
      int di = hs.mesh.extremal_dual(pp);
      auto dp = hs.mesh.dvert[di].head(3) / hs.mesh.dvert[di][3];
      auto d0 = pp.dot(dp);
      auto [pj, v] = hs.mesh.closest_primal(pp);

      mx = std::min(mx, -d0 - os);
      mn = std::min(mn, -v - os);
    }
    if (randomize) {
      double w = dis(gen);
      // res[n] = (w * mn + (1 - w) * mx);
      res[n] = (w * mn + (1 - w) * std::max(mx, mn - top_cell_shift));
    } else {
      // res[n] = 0.5 * (mn + mx);
      res[n] = 0.5 * (mn + std::max(mx, mn - top_cell_shift));
    }
  }

#ifdef PROFILING
  time->end_section();
#endif

  return res;
}


template <typename T>
rasterizer laguerre_diagram<T>::get_rasterizer(std::optional<std::function<Eigen::Vector2d(Eigen::Vector2d)>> transform) {

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
      sidx[i + 1] = sidx[i] + e.ls.size() - 1;
  }

  // rasterizer arguments
  std::vector<std::tuple<Eigen::Vector2d, Eigen::Vector2d, int>> seg;
  auto blo = (transform ? transform.value()(sim.spmin) : sim.spmin);
  auto bhi = (transform ? transform.value()(sim.spmax) : sim.spmax);
  Eigen::Array4d bounds = {
    std::min(blo(0), bhi(0)), std::max(blo(0), bhi(0)),
    std::min(blo(1), bhi(1)), std::max(blo(1), bhi(1))
  };

  bool flipy = (blo(1) > bhi(1));

  // collect segments
  for (int i = 0; i < (int)edglist.size(); ++i) {
    const auto &e = edglist[i];

    if (e.pj == 2 || e.pj == 3)
      continue;

    // get (primal) index of cell below:
    // if the edge is oriented left-to-right then pi is below, if it goes right-to-left then pj is below
    // indices >= n (i.e. the top cell) are all combined into n
    int vidx = std::min((((e.ls.x(0,0) > e.ls.x(e.ls.size()-1,0)) ^ flipy) ? e.pj : e.pi) - 6, n);
    if (vidx < 0) vidx = std::min(((e.ls.x(0,0) > e.ls.x(e.ls.size()-1,0)) ? e.pi : e.pj) - 6, n);
    if (vidx < 0) throw std::runtime_error("in preparing rasterizer, edge has both primal vertices outside: " + e.repr());

    // add segments
    int j = sidx[i];
    for (int i = 1; i < e.ls.size(); ++i, ++j) {
      if ((int)seg.size() != j)
        throw std::runtime_error(FORMAT("segment index does not match expected: {} instead of {}", seg.size(), j));
      
      seg.push_back({
          (transform ? transform.value()(e.ls.x.row(i-1)) : e.ls.x.row(i-1)),
          (transform ? transform.value()(e.ls.x.row(i)) : e.ls.x.row(i)),
          vidx
        });
    }
  }

  // add a segment across the top (for boundary cell)
  seg.push_back({
      transform ? transform.value()(Eigen::Vector2d(sim.spmin(0), sim.spmax(1) + 1)) : Eigen::Vector2d(sim.spmin(0), sim.spmax(1) + 1),
      transform ? transform.value()(Eigen::Vector2d(sim.spmax(0), sim.spmax(1) + 1)) : Eigen::Vector2d(sim.spmax(0), sim.spmax(1) + 1),
      n});

  // add a segment below the bottom
  seg.push_back({
      transform ? transform.value()(Eigen::Vector2d(sim.spmin(0), sim.spmin(1) * 0.5)) : Eigen::Vector2d(sim.spmin(0), sim.spmin(1) * 0.5),
      transform ? transform.value()(Eigen::Vector2d(sim.spmax(0), sim.spmin(1) * 0.5)) : Eigen::Vector2d(sim.spmax(0), sim.spmin(1) * 0.5),
      n});
  
  rasterizer rast(seg, bounds);

#ifdef PROFILING
  time->end_section();
#endif

  return rast;
}

template <typename T>
std::unordered_map<uint64_t, double> laguerre_diagram<T>::jac() {
#ifdef PROFILING
  time->start_section(idx_time_jacobian);
#endif
  // get jacobian (derivative of masses w.r.t. dual)

  std::unordered_map<uint64_t, double> res;
  for (const auto &e : edglist) {
    if (e.pj >= 6) {
      uint64_t i = e.pi-6;
      uint64_t j = std::min(e.pj-6,n);

      res[(i << 32) | i] += e.dif;
      res[(i << 32) | j] -= e.dif;

      res[(j << 32) | j] += e.dif;
      res[(j << 32) | i] -= e.dif;
    }
  }

  for (const auto &e : vedglist) {
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

template <typename T>
std::pair<Eigen::VectorXd, std::pair<Eigen::VectorXi, Eigen::VectorXi>> laguerre_diagram<T>::jac_coo() {
  auto J = jac();

#ifdef PROFILING
  time->start_section(idx_time_jacobian);
#endif

  Eigen::VectorXd data(J.size());
  Eigen::VectorXi i(J.size()), j(J.size());

  int k = 0;
  for (const auto &v : J) {
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

template <typename T>
Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> laguerre_diagram<T>::get_poly(const int &pi_t) {
  int pi = pi_t + 6;
  std::vector<diagram_edge*> edg, path;
  int len = 0;
  for (const auto &ei : faces[pi]) {
    edg.push_back(&edglist[ei]);
    len += edglist[ei].ls.size() - 1;
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
  for (const auto &ep : path)
    if (ep->pi == pi)
      for (int i = 0; i < ep->ls.size() - 1; ++i)
        res.row(ri++) = ep->ls.x.row(i);
    else
      for (int i = ep->ls.size() - 1; i > 0; --i)
        res.row(ri++) = ep->ls.x.row(i);

  return res;
}

template <typename T>
void laguerre_diagram<T>::bind(py::module_ &m) {
  py::class_<laguerre_diagram<T>, std::shared_ptr<laguerre_diagram<T>>>(m, ("LaguerreDiagram_" + type_name<T>::value()).c_str())
    .def(py::init<const laguerre_diagram<T>::seeds_t&,
         const Eigen::Ref<const laguerre_diagram<T>::VecX>&,
         std::shared_ptr<physical_parameters>,
         const simulation_parameters&>(),
         py::arg("ys"), py::arg("duals"),
         py::arg("phys"), py::arg("sim"))
    .def(py::init<std::shared_ptr<laguerre_diagram<T>>,
         const Eigen::Ref<const laguerre_diagram<T>::VecX>&,
         const bool&>(),
         py::arg("parent"), py::arg("duals"), py::arg("keep_parent")=true)
    .def("detach", &laguerre_diagram<T>::detach)
    .def_readwrite("parent", &laguerre_diagram<T>::parent)
    .def_readonly("n", &laguerre_diagram<T>::n)
    .def_readonly("phys", &laguerre_diagram<T>::phys)
    .def_readonly("sim", &laguerre_diagram<T>::sim)
    .def_readonly("hints", &laguerre_diagram<T>::hints)
    .def_readonly("ys", &laguerre_diagram<T>::ys)
    .def_readonly("duals", &laguerre_diagram<T>::duals)
    .def_readonly("hs", &laguerre_diagram<T>::hs)
    .def_readonly("edglist", &laguerre_diagram<T>::edglist)
    .def_readonly("verts", &laguerre_diagram<T>::verts)
    .def_readonly("faces", &laguerre_diagram<T>::faces)
    .def_readonly("dfacets", &laguerre_diagram<T>::dfacets)
    .def_readonly("areas", &laguerre_diagram<T>::areas)
    .def_readonly("areaerrs", &laguerre_diagram<T>::areaerrs)
    .def_readonly("boundary_spt", &laguerre_diagram<T>::boundary_spt)
#ifdef PROFILING
    .def_readonly("time", &laguerre_diagram<T>::time)
#endif
    .def("touching_dual", &laguerre_diagram<T>::touching_dual,
         py::arg("randomize")=false)
    .def("jac_coo", &laguerre_diagram<T>::jac_coo)
    .def("get_rasterizer", &laguerre_diagram<T>::get_rasterizer, py::arg("transform") = std::nullopt)
    .def("get_poly", &laguerre_diagram<T>::get_poly);

  m.def("LaguerreDiagram",
        [](const laguerre_diagram<T>::seeds_t &ys,
           const Eigen::Ref<const laguerre_diagram<T>::VecX> &duals,
           std::shared_ptr<physical_parameters> phys,
           const simulation_parameters &sim) {
          return laguerre_diagram<T>(ys, duals, phys, sim);
        },
        py::arg("ys"), py::arg("duals"),          
        py::arg("phys"), py::arg("sim"));

  m.def("LaguerreDiagram",
        [](std::shared_ptr<laguerre_diagram<T>> parent,
           const Eigen::Ref<const laguerre_diagram<T>::VecX> &duals,
           const bool &keep_parent) {
          return laguerre_diagram<T>(parent, duals, keep_parent);
        },
        py::arg("parent"), py::arg("duals"), py::arg("keep_parent")=true);
  
}
