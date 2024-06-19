// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_SECTIONTREE_H
#define HLBGMYTOOL_OCT_SECTIONTREE_H

#include <fstream>
#include <list>
#include <vector>

#include "FluidSiteTree.h"
#include "MaskTree.h"
#include "Oct.h"

// This is a flattened Octree that at the lowest level stores zero or
// more Sections, in the meaning of PETSc.

// This is a pair of integers (offset and a count) and some data. The
// integers index into the data.

namespace hemelb {

namespace H5 {
class Group;
using GroupPtr = std::shared_ptr<Group>;
}  // namespace H5

namespace gmytool::oct {

template <class T>
struct Section {
  using IndT = uint64_t;

  void append() {
    offsets.push_back(data.size());
    counts.push_back(0);
  }

  template <class... Args>
  void append(Args&&... args) {
    offsets.push_back(data.size());
    counts.push_back(1);
    data.emplace_back(std::forward<Args>(args)...);
  }

  void write(H5::GroupPtr grp) const;

  std::vector<IndT> offsets;
  std::vector<IndT> counts;
  std::vector<T> data;
};

// The tree data is stored in nLevel vectors

// Each branch node in the tree is 8 indices giving the offsets of the
// children in the next level down's vector.

// Leaf nodes are the sections.

// To keep addressing consistent we add an empty vector for the leaf nodes.

class SectionTree {
 public:
  using Ptr = std::shared_ptr<SectionTree>;

  using Int = MaskTree::Int;
  using IndT = uint64_t;

  using TreeLevel = std::vector<IndT>;
  using Tree = std::vector<TreeLevel>;

  static constexpr IndT NA() { return ~0; };

  static inline Int LocalOffset(Int i, Int j, Int k, Int lvl) {
    Int xbit = (i >> lvl) & 1;
    Int ybit = (j >> lvl) & 1;
    Int zbit = (k >> lvl) & 1;
    return (xbit << 2) | (ybit << 1) | zbit;
  }

  IndT FindIndex(Int i, Int j, Int k) const;

  // void AddSection(std::string& name);

  const Tree& GetTree() const;

  void Write(const std::string& fn) const;

 private:
  friend class SectionTreeBuilder;
  friend class SectionTreeTests;

  SectionTree(size_t nl);

  Int nLevels;
  Tree indices;
  Tree counts;
  IndT total;

  Section<SVector> wall_normals;
  Section<std::array<Link, 26>> links;
};

}  // namespace gmytool::oct
}  // namespace hemelb
#endif
