// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <boost/lockfree/queue.hpp>
#include <exception>
#include "range.hpp"

#include "FloodFill.h"

#include "Neighbours.h"
#include "util/Vector3D.h"
#define unpack(idx) idx[0], idx[1], idx[2], idx[3]

namespace hemelb::gmytool::oct {

FloodFill::FloodFill(const FluidTree& t) : tree(t) {}

auto FloodFill::GetStart() const -> Idx {
  Idx ans;
  try {
    tree.IterDepthFirst(0, 0, [&ans](FluidTree::ConstNodePtr n) {
      ans = {n->X(), n->Y(), n->Z(), 0};
      throw StopIteration();
    });
  } catch (StopIteration& stop) {
    return ans;
  }
  throw std::runtime_error("Could not get a start point from tree");
}

#define PrintIdx(i)                                                        \
  std::cout << "[" << i[0] << ", " << i[1] << ", " << i[2] << ", " << i[3] \
            << "]" << std::endl

// Helper functor that will create a leaf node and increment the fluid
// count all the way up to the root. User must be sure that it will
// only be called on coordinates that don't exist else the fluid
// counts will be wrong!!
class CreateAndIncrementor {
 public:
  typedef MaskTree::NodePtr NodePtr;
  typedef MaskTree::Int Int;

  CreateAndIncrementor(MaskTree& t) : tree(t) {}

  void operator()(Int i, Int j, Int k) {
#ifndef NDEBUG
    // Only unused in production mode
    unsigned initial = tree.Root()->Data();
#endif
    auto path = tree.Root()->GetCreatePath(i, j, k, 0);
    for (auto n : path) {
      n->Data() += 1;
    }
    assert(tree.Root()->Data() == (initial + 1));
  }

 private:
  MaskTree& tree;
};

MaskTree FloodFill::operator()() const {
  auto constexpr& dirs = Neighbours::Displacements;
  typedef boost::lockfree::queue<Idx> Queue;
  auto seed = GetStart();

  MaskTree seen(tree.Level());
  // seen.Root()->Data() = 0;
  CreateAndIncrementor candi(seen);
  candi(seed[0], seed[1], seed[2]);

  // WorkQ holds sites that are def fluid but we don't know about their
  // neighbours
  Queue workQ(100);
  workQ.unsynchronized_push(seed);

  auto BULK = std::make_shared<FluidSite>();

  while (!workQ.empty()) {
    // Grab the point from the queue
    Idx pt;
    if (!workQ.pop(pt))
      continue;

    auto leaf_node_ptr = tree.Get(unpack(pt));

    // enqueue (pt + dirs[i]) if we have not seen it
    auto enqueue_neighbour_if_unseen = [&](unsigned i) {
      Idx neigh{FluidTree::Int(pt[0] + dirs[i][0]),
                FluidTree::Int(pt[1] + dirs[i][1]),
                FluidTree::Int(pt[2] + dirs[i][2]), pt[3]};

      if (!seen.Get(unpack(neigh))) {
        candi(neigh[0], neigh[1], neigh[2]);
        workQ.push(neigh);
      }
    };

    if (leaf_node_ptr) {
      // Node exists already: it must be fluid and have been created
      // by SurfaceVoxeliser. Look at all links and queue all sites
      // with intersection == none
      const auto& node_data = leaf_node_ptr->Data();
      const auto& links = node_data.leaf->links;
      for (auto i : range(0, 26)) {
        if (links[i].type == Intersection::None) {
          enqueue_neighbour_if_unseen(i);
        }
      }
    } else {
      // Node does not exist, but it must be fluid (as we are only
      // enqueuing nodes that are linked to a fluid site without cuts)

      // Could create a leaf node for it with the basic fluid site
      // state, but that is a lot of work

      // Instead, just iter over all the neighbours and enqueue if not
      // seen
      for (auto i : range(0, 26)) {
        enqueue_neighbour_if_unseen(i);
      }
    }
  }
  return seen;
}
}  // namespace hemelb::gmytool::oct
