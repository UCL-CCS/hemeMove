// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "SectionTreeBuilder.h"

namespace hemelb::gmytool::oct {

SectionTreeBuilder::SectionTreeBuilder(const MaskTree& mask,
                                       const FluidTree& edges)
    : maskTree(mask),
      edgeTree(edges),
      nLevels(mask.Level()),
      nEdgeSites(edges.Root()->Data().count),
      offsets(mask.Level() + 1),
      edge_ptrs(mask.Level() + 1),
      wall_normal_offset(0),
      links_offset(0) {}

SectionTree::Ptr SectionTreeBuilder::operator()() {
  output = SectionTree::Ptr(new SectionTree(nLevels));
  maskTree.Root()->Accept(*this);
  return output;
}

// edge_ptrs is a vector of tree nodes with the path to the current
// node in the edge tree (or null if no corresponding one)
void SectionTreeBuilder::Arrive(MaskTree::ConstNodePtr np) {
  const MaskTree::Node& n = *np;
  Int lvl = n.Level();

  // First, get our equivalent in the edge tree
  if (lvl == nLevels) {
    // We are the root node - easy!
    edge_ptrs[lvl] = edgeTree.Root();
  } else {
    // branch or leaf - must check parent equivalent is present as could be null
    auto edge_parent = edge_ptrs[lvl + 1];
    if (edge_parent) {
      // It exists so try to find our equivalent
      edge_ptrs[lvl] = edge_parent->Get(n.X(), n.Y(), n.Z(), lvl);
    } else {
      // Parent doesn't exist so our equiv doesnt either
      edge_ptrs[lvl] = nullptr;
    }
  }

  // This is the first we knew of this node's existence (unless we are
  // root), so must store its position into the parent's data
  if (lvl < nLevels) {
    auto parent_current_i = offsets[lvl + 1];
    auto loffset = LocalOffset(n);

    // Write this level's offset value into parent's current data oct
    // We will then increment offset[lvl] on Departure.
    output->indices[lvl + 1][parent_current_i + loffset] = offsets[lvl];

  } else {
    // nothing
  }

  if (lvl) {
    // We are not a leaf node so we must allocate space for our
    // potential children
    auto new_size = output->indices[lvl].size() + 8;

    // Fill with null children until they add themselves
    output->indices[lvl].resize(new_size, SectionTree::NA());
    // Or zeros in case of counts
    output->counts[lvl].resize(new_size, 0);
  } else {
    // We are at a leaf node - we must extend the sections
    if (edge_ptrs[0]) {
      // We are an edge site so we have data to add
      auto& edge_site = *(edge_ptrs[0]->Data().leaf);

      output->links.append(edge_site.links);
      if (edge_site.has_normal) {
        output->wall_normals.append(edge_site.normal);
      } else {
        output->wall_normals.append();
      }

    } else {
      // we are a bulk site
      output->links.append();
      output->wall_normals.append();
    }
  }
}

void SectionTreeBuilder::Depart(MaskTree::ConstNodePtr n) {
  Int lvl = n->Level();

  // Store contained site count
  SectionTree::IndT site_count = 0;
  if (lvl) {
    // non-leaf so sum my child counts
    auto iter = output->counts[lvl].cbegin() + offsets[lvl];
    for (auto i = 0; i < 8; ++i) {
      site_count += *iter;
    }
  } else {
    site_count = 1;
  }

  if (lvl < nLevels) {
    // Non-root
    auto loffset = LocalOffset(*n);
    output->counts[lvl + 1][offsets[lvl + 1] + loffset] = site_count;
  } else {
    // Root
    output->total = site_count;
  }

  // Update offset (1 for leaf nodes)
  offsets[lvl] += lvl ? 8 : 1;

  // Clear the edge pointer to make sure we don't leave it hanging around
  edge_ptrs[lvl] = nullptr;
}

}  // namespace hemelb::gmytool::oct
