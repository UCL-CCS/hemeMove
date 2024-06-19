// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "InconsistentFluidnessError.h"
#include <sstream>
#include "Site.h"

namespace hemelb::gmytool::gmy {

namespace {

void FormatSite(std::ostringstream& msg, const Site& site) {
  msg << "site index " << site.GetIndex() << ", position " << site.Position
      << ", which is ";
  if (site.IsFluid)
    msg << "fluid";
  else
    msg << "solid";
}

void FormatBoolInside(std::ostringstream& msg, const bool sinside) {
  if (sinside == true)
    msg << "fluid";
  else
    msg << "solid";
}

}  // namespace

InconsistentFluidnessError::InconsistentFluidnessError(const Site& s1,
                                                       const Site& s2,
                                                       const int nHits)
    : s1(s1), s2(s2), nHits(nHits) {
  std::ostringstream msgStream;
  msgStream << "Inconsistent fluidness detected between ";
  FormatSite(msgStream, s1);
  msgStream << " and ";
  FormatSite(msgStream, s2);
  msgStream << " but found " << nHits << " intersections with the surface.";
  this->msg = msgStream.str();
}

InconsistentIntersectRayError::InconsistentIntersectRayError(const Site& s1,
                                                             const Site& s2,
                                                             const int nHits,
                                                             const bool sinside,
                                                             const bool ninside)
    : s1(s1), s2(s2), nHits(nHits), sinside(sinside), ninside(ninside) {
  std::ostringstream msgStream;
  msgStream << "Inconsistent fluidness from ray and intersect between ";
  FormatSite(msgStream, s1);
  msgStream << " and ";
  FormatSite(msgStream, s2);
  msgStream
      << " found " << nHits
      << " intersections with the surface, but according to Raycast site 1 is ";
  FormatBoolInside(msgStream, sinside);
  msgStream << " while site 2 is ";
  FormatBoolInside(msgStream, ninside);
  this->msg = msgStream.str();
}

const char* InconsistentFluidnessError::what() const noexcept {
  return this->msg.c_str();
}

const char* InconsistentIntersectRayError::what() const noexcept {
  return this->msg.c_str();
}
}  // namespace hemelb::gmytool::gmy
