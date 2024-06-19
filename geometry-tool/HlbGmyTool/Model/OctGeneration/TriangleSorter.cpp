// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "TriangleSorter.h"
#include <algorithm>
#include <deque>
#include <future>

namespace hemelb::gmytool::oct {

template <>
void WriteNodeData<IdList>(std::ostream& os, const IdList& data) {
  os << '[';
  bool first = true;
  for (auto id : data) {
    if (!first)
      os << ", ";

    os << id;
    first = false;
  }
  os << ']';
}

TriTree TrianglesToTreeSerial(const int n_levels,
                              const int tri_level,
                              const std::vector<Vector>& points,
                              const std::vector<Index>& triangles) {
  return TrianglesToTree_Worker(n_levels, tri_level, points, triangles.begin(),
                                triangles.end(), 0);
}

TriTree TrianglesToTree_Worker(const int n_levels,
                               const int tri_level,
                               const std::vector<Vector>& points,
                               std::vector<Index>::const_iterator triStart,
                               std::vector<Index>::const_iterator triEnd,
                               const int tri_index_start) {
  typedef TriTree::Int Int;
  Int cube_size = 1 << n_levels;
  Int tri_box_size = 1 << tri_level;

  std::vector<Int> non_halo_edges, lowers, uppers;
  for (Int i = 0; i < cube_size; i += tri_box_size) {
    non_halo_edges.push_back(i);
    // i==0 use lower=0 cos we know that no tri can be outside our bounding box.
    // The upper and lower bounds of the regions of influence of the tri_level
    // nodes.
    // +/- 1 because we have to consider the halo.
    lowers.push_back(i ? i - 1 : 0);
    uppers.push_back(i + tri_box_size);
  }

  TriTree tree(n_levels);
  // for (auto i_tri = 0; i_tri < ntris; ++i_tri) {
  for (auto tri_it = triStart; tri_it != triEnd; ++tri_it) {
    auto triPtIds = *tri_it;

    // bounding box of the triangle
    Vector min(std::numeric_limits<double>::infinity());
    Vector max(0);
    for (auto iPt : triPtIds) {
      min.UpdatePointwiseMin(points[iPt]);
      max.UpdatePointwiseMax(points[iPt]);
    }

    Index range_min;
    Index range_max;
    for (auto iDim = 0; iDim < 3; ++iDim) {
      // Figure out the first tri box to add this tri to
      // iterator pointing to the first element > min
      auto upper_ptr =
          std::upper_bound(uppers.cbegin(), uppers.cend(), min[iDim]);
      // iterator pointing to the first element > max
      auto lower_ptr =
          std::upper_bound(lowers.cbegin(), lowers.cend(), max[iDim]);
      // convert to indices
      range_min[iDim] = upper_ptr - uppers.cbegin();
      range_max[iDim] = lower_ptr - lowers.cbegin();
    }

    // This bit is the only non-thread-safe part.
    // We will later merge the trees if this was run in parallel
    for (auto i = range_min[0]; i < range_max[0]; ++i)
      for (auto j = range_min[1]; j < range_max[1]; ++j)
        for (auto k = range_min[2]; k < range_max[2]; ++k) {
          auto node = tree.GetCreate(non_halo_edges[i], non_halo_edges[j],
                                     non_halo_edges[k], tri_level);
          node->Data().insert(tri_it - triStart + tri_index_start);
        }
  }
  return tree;
}

TreeSummer::TreeSummer(TriTree::Int l, TriTree::Int tl)
    : levels(l), tri_level(tl), tree(l) {}

TriTree& TreeSummer::GetTree() {
  return tree;
}

void TreeSummer::Add(TriTree& source) {
  source.IterDepthFirst(tri_level, tri_level, [&](TriTree::NodePtr src) {
    auto dest = tree.GetCreate(src->X(), src->Y(), src->Z(), tri_level);
    dest->Data().insert(boost::container::ordered_unique_range_t(),
                        src->Data().begin(), src->Data().end());
  });
}

TriTree TrianglesToTreeParallel(const int n_levels,
                                const int tri_level,
                                const std::vector<Vector>& points,
                                const std::vector<Index>& triangles,
                                int nprocs) {
  if (nprocs == 0)
    nprocs = std::thread::hardware_concurrency();

  auto n_tri = triangles.size();
  // tris per proc
  auto tpp = float(n_tri) / float(nprocs);

  // procs = []
  // # answer queue
  // q = multiprocessing.Queue()
  std::deque<std::future<TriTree>> futures;
  for (auto i = 0; i < nprocs; ++i) {
    decltype(n_tri) min_ind = i * tpp;
    decltype(n_tri) max_ind = std::min(decltype(n_tri)((i + 1) * tpp), n_tri);
    auto triStart = triangles.begin() + min_ind;
    auto triEnd = triangles.begin() + max_ind;
    futures.push_back(std::async(std::launch::async, TrianglesToTree_Worker,
                                 n_levels, tri_level, points, triStart, triEnd,
                                 min_ind));
  }

  TreeSummer summer(n_levels, tri_level);

  while (!futures.empty()) {
    // Future's are non-copyable - have to move them in & out of the queue
    auto fut = std::move(futures.front());
    futures.pop_front();
    auto status = fut.wait_for(std::chrono::milliseconds(1));
    if (status == std::future_status::ready) {
      auto tree = fut.get();
      summer.Add(tree);
    } else {
      futures.push_back(std::move(fut));
    }
  }
  return summer.GetTree();
}

}  // namespace hemelb::gmytool::oct
