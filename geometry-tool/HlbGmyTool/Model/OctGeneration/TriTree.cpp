// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "TriTree.h"

namespace hemelb::gmytool::oct {

std::ostream& operator<<(std::ostream& os, const IdList& lst) {
  bool first = true;
  os << '[';
  for (auto i : lst) {
    if (!first)
      os << ',';
    os << i;
    first = false;
  }
  os << ']';
  return os;
}

}  // namespace hemelb::gmytool::oct
