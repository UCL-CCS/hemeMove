// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_NEIGHBOURS_H
#define HLBGMYTOOL_GMY_NEIGHBOURS_H

#include <vector>
#include "Index.h"

#include "io/formats/geometry.h"

namespace hemelb::gmytool::gmy {

struct Neighbours {
  // shortcut to geometry class
  using gmy = hemelb::io::formats::geometry;

  // This assumes that the hemelb lattice descriptor class has a zero vector.
  static constexpr int n = gmy::NumberOfDisplacements;
  static constexpr int nLater = n / 2;

  static std::vector<unsigned int> laterNeighbourIndices;
  static std::vector<gmy::Displacement> vectors;
  static std::vector<double> norms;
  static std::vector<unsigned int> inverses;

  static void Init();

 private:
  // Private to ensure it's not instantiated.
  Neighbours();
};

}  // namespace hemelb::gmytool::gmy
#endif  // HLBGMYTOOL_GMY_NEIGHBOURS_H
