// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_TRIANGLESORTER_H
#define HLBGMYTOOL_OCT_TRIANGLESORTER_H

#include <vector>
#include "Index.h"
#include "TriTree.h"

namespace hemelb::gmytool::oct {

// Create an Octree with n_levels potential number of levels.
// Nodes will be created down to the level tri_level but only if voxel-voxel
// links originating within that node could potentially intersect the triangles
// specified by (points, triangles), where:
// points is a numpy array of doubles with shape == (nPoints, 3) specifying the
// vertices. triangles is a numpy array of ints with shape == (nTris, 3)
// specifiying the IDs of the verticies comprising the triangle. The create
// nodes at tri_level will have an attribute 'triIds' added which is a list of
// the IDs of the triangles that might intersect its voxel-voxel links.

TriTree TrianglesToTreeSerial(const int n_levels,
                              const int tri_level,
                              const std::vector<Vector>& points,
                              const std::vector<Index>& triangles);

TriTree TrianglesToTreeParallel(const int n_levels,
                                const int tri_level,
                                const std::vector<Vector>& points,
                                const std::vector<Index>& triangles,
                                int nprocs);

TriTree TrianglesToTree_Worker(const int n_levels,
                               const int tri_level,
                               const std::vector<Vector>& points,
                               std::vector<Index>::const_iterator triStart,
                               std::vector<Index>::const_iterator triEnd,
                               const int tri_index_start);

// This guy accumulates multiple trees together, merging their id lists
class TreeSummer {
 public:
  TreeSummer(TriTree::Int levels, TriTree::Int tri_level);
  void Add(TriTree& source);
  TriTree& GetTree();

 private:
  const TriTree::Int levels;
  const TriTree::Int tri_level;
  // This will be the output
  TriTree tree;
};

}  // namespace hemelb::gmytool::oct
#endif
