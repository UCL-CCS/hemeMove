// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_SURFACEVOXELISER_H
#define HLBGMYTOOL_OCT_SURFACEVOXELISER_H

#include <array>
#include <memory>

#include "Cgal.h"
#include "CgalSearch.h"
#include "FluidSiteTree.h"
#include "Index.h"
#include "Iolet.h"
#include "Oct.h"
#include "TriTree.h"

namespace hemelb::gmytool::oct {

struct Cut {
  Cut() : dist(std::numeric_limits<double>::infinity()), id(-1) {}
  double dist;
  int id;
};

// This represents a site (whether fluid or solid) that has one or more
// links cutting the surface.
struct EdgeSite {
  std::array<Cut, 26> closest_cut;
};

class SurfaceVoxeliser {
 public:
  typedef FluidTree::Int Int;

  SurfaceVoxeliser(const int node_size,
                   const std::vector<Vector>& points,
                   const std::vector<Index>& triangles,
                   const std::vector<Vector>& normals,
                   const std::vector<int>& labels,
                   const std::vector<Iolet>& iolets);

  // We maybe need a destructor because the order of destruction of mesh and
  // searcher seems to matter?
  //~SurfaceVoxeliser();

  void ComputeIntersectionsForSite(Int x,
                                   Int y,
                                   Int z,
                                   EdgeSite& outleaf) const;
  std::unique_ptr<FluidSite> ClassifySite(const EdgeSite& node) const;

  FluidTree operator()(const TriTree& inTree, const TriTree::Int tri_level);

 private:
  const static int NDIR = 13;
  const std::vector<Vector>& Points;
  const std::vector<Index>& Triangles;
  const std::vector<Vector>& Normals;
  const std::vector<int>& Labels;
  const std::vector<Iolet>& Iolets;

  CgalMeshPtr mesh;
  std::unique_ptr<CgalSearchTree> searcher;
};

}  // namespace hemelb::gmytool::oct
#endif
