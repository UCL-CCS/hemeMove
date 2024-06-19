// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_TRITREE_H
#define HLBGMYTOOL_OCT_TRITREE_H

#include <boost/container/flat_set.hpp>
#include "Oct.h"

namespace hemelb::gmytool::oct {

typedef boost::container::flat_set<int> IdList;
typedef Octree<IdList> TriTree;

std::ostream& operator<<(std::ostream& os, const IdList& lst);
}  // namespace hemelb::gmytool::oct

#endif
