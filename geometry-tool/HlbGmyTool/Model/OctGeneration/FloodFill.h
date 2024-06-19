// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_FLOODFILL_H
#define HLBGMYTOOL_OCT_FLOODFILL_H

#include "FluidSiteTree.h"
#include "MaskTree.h"

namespace hemelb::gmytool::oct {

// This function flood fills the domain
// And (eventually) propagates the fluid site count of children up the
// tree
class FloodFill {
 public:
  typedef std::array<FluidTree::Int, 4> Idx;

  FloodFill(const FluidTree& t);

  MaskTree operator()() const;

  Idx GetStart() const;

 private:
  class StopIteration : public std::exception {};

  const FluidTree& tree;
};

}  // namespace hemelb::gmytool::oct
#endif
