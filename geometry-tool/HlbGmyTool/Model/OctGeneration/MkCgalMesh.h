// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_MKCGALMESH_H
#define HLBGMYTOOL_OCT_MKCGALMESH_H

#include <vector>
#include "Cgal.h"
#include "Index.h"

namespace hemelb::gmytool::oct {
// Create a CGAL 3D Polyhedral surface based on the input
// Labels each face with its index in the input - use that to index the
// surface type vector
CgalMeshPtr MkCgalMesh(const std::vector<Vector>& ptsIn,
                       const std::vector<Index>& polysIn);
}  // namespace hemelb::gmytool::oct
#endif  // HLBGMYTOOL_OCT_MKCGALMESH_H
