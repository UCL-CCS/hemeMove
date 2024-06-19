// -*- mode: C++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_TESTRESOURCES_HELPERS_H
#define HLBGMYTOOL_OCT_TESTRESOURCES_HELPERS_H

namespace hemelb::gmytool::oct {

template <class... Types>
std::string GetResource(const Types&... paths);

// Zero arg case gets the path to folder from environment.
template <>
inline std::string GetResource() {
  if (const char* env_dir = std::getenv("TESTRESOURCES")) {
    return std::string{env_dir};
  } else {
    throw std::runtime_error("Must set TESTRESOURCES env var");
  }
}

template <class T>
std::string GetResource(const T& pth) {
  return GetResource() + "/" + pth;
}

template <class T, class... Types>
std::string GetResource(const T& pth, const Types&... paths) {
  return GetResource(pth) + "/" + GetResource(paths...);
}

}  // namespace hemelb::gmytool::oct
#endif
