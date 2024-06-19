// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "Neighbours.h"

#include "io/formats/geometry.h"
// shortcut to geometry class

namespace hemelb::gmytool::gmy {

void Neighbours::Init() {
  // copy in the vectors of the neighbourhood
  Neighbours::vectors = {gmy::Neighbourhood.begin(), gmy::Neighbourhood.end()};

  // unsigned int iLater = 0;
  for (unsigned int i = 0, iLater = 0; i < Neighbours::n; ++i) {
    Neighbours::norms[i] =
        std::sqrt(Neighbours::vectors[i].GetMagnitudeSquared());

    // Figure out if this takes us to a later site
    if (Neighbours::vectors[i][0] > 0) {
      // yes
      Neighbours::laterNeighbourIndices[iLater] = i;
      iLater++;
    } else if (Neighbours::vectors[i][0] == 0) {
      // Maybe
      if (Neighbours::vectors[i][1] > 0) {
        // yes
        Neighbours::laterNeighbourIndices[iLater] = i;
        iLater++;
      } else if (Neighbours::vectors[i][1] == 0) {
        // Maybe
        if (Neighbours::vectors[i][2] > 0) {
          // yes
          Neighbours::laterNeighbourIndices[iLater] = i;
          iLater++;
        } else {
          // No
        }
      } else {
        // No
      }
    } else {
      // No
    }

    // Complete the table of inverse indices.
    for (unsigned int j = i; j < Neighbours::n; ++j) {
      if (Neighbours::vectors[i] == -Neighbours::vectors[j]) {
        Neighbours::inverses[i] = j;
        Neighbours::inverses[j] = i;
        break;
      }
    }
  }
}

std::vector<unsigned int> Neighbours::laterNeighbourIndices(Neighbours::nLater);
std::vector<Index> Neighbours::vectors(Neighbours::n);
std::vector<unsigned int> Neighbours::inverses(Neighbours::n);
std::vector<double> Neighbours::norms(Neighbours::n);

}  // namespace hemelb::gmytool::gmy
