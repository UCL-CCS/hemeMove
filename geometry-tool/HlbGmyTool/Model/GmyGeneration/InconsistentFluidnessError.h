// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_INCONSISTENTFLUIDNESSERROR_H
#define HLBGMYTOOL_GMY_INCONSISTENTFLUIDNESSERROR_H

#include <string>
#include "GenerationError.h"

namespace hemelb::gmytool::gmy {

class Site;

struct InconsistentFluidnessError : public GenerationError {
  InconsistentFluidnessError(const Site& s1, const Site& s2, const int nHits);
  ~InconsistentFluidnessError() {}
  virtual const char* what() const noexcept;
  const Site& s1;
  const Site& s2;
  const int nHits;

 private:
  std::string msg;
};

struct InconsistentIntersectRayError : public GenerationError {
  InconsistentIntersectRayError(const Site& s1,
                                const Site& s2,
                                const int nHits,
                                const bool sinside,
                                const bool ninside);
  ~InconsistentIntersectRayError() {}
  virtual const char* what() const noexcept;
  const Site& s1;
  const Site& s2;
  const int nHits;
  const bool sinside;
  const bool ninside;

 private:
  std::string msg;
};

}  // namespace hemelb::gmytool::gmy

#endif  // HLBGMYTOOL_GMY_INCONSISTENTFLUIDNESSERROR_H
