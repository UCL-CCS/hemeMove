// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_IOLET_H
#define HLBGMYTOOL_GMY_IOLET_H

#include "Index.h"

namespace hemelb::gmytool::gmy {

struct Iolet {
  Vector Centre;
  Vector Normal;

  double Radius;
  int Id;
  bool IsInlet;
};

}  // namespace hemelb::gmytool::gmy
#endif  // HLBGMYTOOL_GMY_IOLET_H
