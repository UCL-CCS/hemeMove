// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <cassert>
#include <deque>
#include <memory>
#include "range.hpp"

#include "MkCgalMesh.h"
#include "Neighbours.h"
#include "ParallelApply.h"
#include "SegmentFactory.h"
#include "SurfaceVoxeliser.h"

#include <boost/function_output_iterator.hpp>

namespace hemelb::gmytool::oct {

SurfaceVoxeliser::SurfaceVoxeliser(const int ns,
                                   const std::vector<Vector>& p,
                                   const std::vector<Index>& t,
                                   const std::vector<Vector>& n,
                                   const std::vector<int>& l,
                                   const std::vector<Iolet>& io)
    : Points(p),
      Triangles(t),
      Normals(n),
      Labels(l),
      Iolets(io),
      mesh(MkCgalMesh(p, t)),
      searcher(
          new CgalSearchTree(faces(*mesh).first, faces(*mesh).second, *mesh)) {}

using Segment_intersection = boost::optional<
    CgalSearchTree::Intersection_and_primitive_id<CgalSegment>::Type>;
using dist_isec = std::pair<double, int>;

// Helper class for filtering out clusters of intersections
struct Clusterer {
  const std::vector<Vector>& normals;
  std::vector<dist_isec>& output;
  const Vector& direction;

  int n_pts;
  int n_pos;
  int n_neg;
  const dist_isec* current;

  Clusterer(const std::vector<Vector>& n,
            std::vector<dist_isec>& o,
            const Vector& d)
      : normals(n),
        output(o),
        direction(d),
        n_pts(0),
        n_pos(0),
        n_neg(0),
        current(nullptr) {}

  const double& GetEnd() { return current->first; }

  void Reset() {
    n_pts = 0;
    n_pos = 0;
    n_neg = 0;
    current = nullptr;
  }

  void Start(const dist_isec& di) {
    Reset();
    Extend(di);
    current = &di;
  }

  void Extend(const dist_isec& di) {
    n_pts++;
    auto dp = Vector::Dot(normals[di.second], direction);
    if (dp > 0.0) {
      n_pos++;
    } else {
      n_neg++;
    }
    current = &di;
  }

  void Finish() {
    if (n_pos && n_neg) {
      // do not add anything if
    } else {
      output.push_back(*current);
    }
  }
};

// This filters a list of intersections for clusters of intersections that
// should be merged
std::vector<dist_isec> FilterClusters(const CgalPoint& start,
                                      const Vector direction,
                                      const std::vector<Vector>& normals,
                                      std::vector<dist_isec>& dirty,
                                      double tol = 1e-3) {
  std::vector<dist_isec> clean;
  if (dirty.size() > 1) {
    // Look for "clusters" of points where they are really close together

    // Cluster intersections should either be ignored (if a mix of
    // inner and outer) or merged to a single intersection (with an
    // arbitrary face picked).
    bool first = true;
    clean.reserve(dirty.size());
    Clusterer clst(normals, clean, direction);
    for (const auto& cur : dirty) {
      if (first) {
        clst.Start(cur);
        first = false;
        continue;
      }

      if (cur.first < (clst.GetEnd() + tol)) {
        clst.Extend(cur);
      } else {
        clst.Finish();
        clst.Start(cur);
      }
    }
    clst.Finish();

  } else {
    // zero or one intersections: just swap vectors
    std::swap(clean, dirty);
  }
  return clean;
}

void SurfaceVoxeliser::ComputeIntersectionsForSite(Int x,
                                                   Int y,
                                                   Int z,
                                                   EdgeSite& outleaf) const {
  auto constexpr& directions = Neighbours::Displacements;
  Index siteIdx(x, y, z);

  for (unsigned iDir = 0; iDir < directions.size(); ++iDir) {
    auto const& direction = directions[iDir];
    double inorm = 1.0 / std::sqrt(direction.GetMagnitudeSquared());
    auto neighIdx = siteIdx + direction;

    CgalPoint start(siteIdx.x, siteIdx.y, siteIdx.z);
    CgalPoint end(neighIdx.x, neighIdx.y, neighIdx.z);
    CgalSegment link(start, end);

    // This will hold all dist, intersection pairs
    std::vector<dist_isec> intersections;
    // Computes dist
    auto dist_from_start = [&start, &inorm](const CgalPoint& pt) -> double {
      return std::sqrt(CGAL::squared_distance(start, pt)) * inorm;
    };
    // Lambda will add intersections keeping them sorted
    auto add_isec = boost::make_function_output_iterator(
        [&intersections, &dist_from_start](const Segment_intersection& isec) {
          try {
            // Calculate the cut distance
            auto d = dist_from_start(boost::get<CgalPoint>(isec->first));
            // Figure out where to insert it
            auto insertion_point = std::lower_bound(
                intersections.begin(), intersections.end(), d,
                [](const dist_isec& d_i, const double& d) -> bool {
                  return d_i.first < d;
                });
            auto& id = isec->second->id();
            intersections.emplace(insertion_point, d, id);
          } catch (boost::bad_get& e) {
            // Segment - ignore
          }
        });
    // Run
    searcher->all_intersections(link, add_isec);

    // Clean up clusters
    auto clean_intersections =
        FilterClusters(start, direction.as<double>(), Normals, intersections);

    // Get the closest and add if it exists
    auto closest = clean_intersections.begin();
    if (closest != clean_intersections.end()) {
      auto& cuts = outleaf.closest_cut;
      cuts[iDir].dist = closest->first;
      cuts[iDir].id = closest->second;
    }
  }
}

std::unique_ptr<FluidSite> SurfaceVoxeliser::ClassifySite(
    const EdgeSite& node) const {
  // Create output subtree
  auto constexpr& dirs = Neighbours::Displacements;

  auto& cuts = node.closest_cut;
  // Examine all the links totting up how many are different types
  unsigned nIn = 0;
  unsigned nOut = 0;
  unsigned nNA = 0;

  auto fsite = std::make_unique<FluidSite>();
  for (auto i : range(26)) {
    if (cuts[i].dist < 1) {
      // there's a cut
      Vector dir{dirs[i]};
      auto& faceId = cuts[i].id;
      auto& norm = this->Normals[faceId];
      auto dp = Vector::Dot(dir, norm);
      if (dp > 0) {
        // We are inside
        ++nIn;

        auto& label = this->Labels[faceId];
        fsite->links[i].dist = cuts[i].dist;
        if (label < 0) {
          // Wall
          fsite->links[i].type = Intersection::Wall;
        } else {
          // Iolet
          // TODO: implement this properly
          auto& iolet = this->Iolets[label];
          fsite->links[i].type =
              iolet.IsInlet ? Intersection::Inlet : Intersection::Outlet;
          fsite->links[i].id = iolet.Id;
        }
      } else {
        ++nOut;
        // This site better not have any inside!
      }
    } else {
      // no cut
      ++nNA;
      fsite->links[i].type = Intersection::None;
    }
  }

  // Now some sanity checks
  assert((nIn + nOut + nNA) == 26);

  if (nIn) {
    // One or more interior sites, this one must be added to output
    if (nOut) {
      // std::cout << node << std::endl;
    }
    assert(nOut == 0);

    // Find closest wall link
    float wall_dist = std::numeric_limits<float>::infinity();
    int faceId;
    for (auto i : range(26)) {
      if (fsite->links[i].type == Intersection::Wall) {
        if (cuts[i].dist < wall_dist) {
          wall_dist = cuts[i].dist;
          faceId = cuts[i].id;
        }
      }
    }
    if (wall_dist < 1) {
      fsite->has_normal = true;
      fsite->normal = SVector{this->Normals[faceId]};
    } else {
      fsite->has_normal = false;
    }

  } else {
    fsite = nullptr;
  }
  return fsite;
}

// This will voxelise the entire input tree.
// Nodes at tri_level will be treated as solid blocks.
// Nodes at that level will be handled in parallel.
FluidTree SurfaceVoxeliser::operator()(const TriTree& inTree,
                                       const TriTree::Int tri_level) {
  using Pool = ParallelApply<FluidTree::NodePtr, TriTree::ConstNodePtr>;
  using Coord = hemelb::util::Vector3D<Int>;

  // Set up pool to call block processing lambda on each tri_level node
  Pool pool(
      [this](TriTree::ConstNodePtr in) {
        Coord start(in->X(), in->Y(), in->Z());
        Int const bsz = 1 << in->Level();
        Coord stop(in->X() + bsz, in->Y() + bsz, in->Z() + bsz);
        auto ans = FluidTree::NodePtr(
            new FluidTree::Branch(in->X(), in->Y(), in->Z(), in->Level()));
        unsigned nleaf = 0;

        // Iter over all coords in input node
        for (auto x : range(start.x, stop.x))
          for (auto y : range(start.y, stop.y))
            for (auto z : range(start.z, stop.z)) {
              EdgeSite intersec_data;
              ComputeIntersectionsForSite(x, y, z, intersec_data);
              auto outleafdata = ClassifySite(intersec_data);
              if (outleafdata) {
                auto outleaf = ans->GetCreate(x, y, z, 0);
                outleaf->Data().leaf = std::move(outleafdata);
                outleaf->Data().count = 1;
                ++nleaf;
              }
            }
        ans->Data().count = nleaf;
        return ans;
      },
      0, std::numeric_limits<int>::max());

  // Submit all the tri_level blocks for processing
  // Merging of output done below
  std::deque<Pool::Future> result_queue;
  inTree.IterDepthFirst(tri_level, tri_level, [&](TriTree::ConstNodePtr node) {
    result_queue.emplace_back(pool(node));
  });
  pool.Done();

  // Merge blocks
  FluidTree outTree(inTree.Level());

  unsigned nLeafNodes = 0;
  while (!result_queue.empty()) {
    auto& fut = result_queue.front();

    auto edge_node = fut.get();
    nLeafNodes += edge_node->Data().count;
    outTree.Set(edge_node->X(), edge_node->Y(), edge_node->Z(), tri_level,
                edge_node);
    result_queue.pop_front();
  }
  outTree.Root()->Data().count = nLeafNodes;
  return outTree;
}

}  // namespace hemelb::gmytool::oct
