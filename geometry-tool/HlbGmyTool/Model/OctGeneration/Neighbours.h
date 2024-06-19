// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_OCT_NEIGHBOURS_H
#define HLBGMYTOOL_OCT_NEIGHBOURS_H

#include "Index.h"
#include "io/formats/geometry.h"
#include "range.hpp"

namespace hemelb::gmytool::oct {

namespace detail {
template <std::size_t N>
static constexpr std::array<Vector, N> i2d(std::array<Index, N> const& idata) {
  std::array<Vector, N> ans;
  for (auto i : range(N)) {
    ans[i].x = idata[i].x;
    ans[i].y = idata[i].y;
    ans[i].z = idata[i].z;
  }
  return ans;
}

template <std::size_t N>
constexpr std::array<unsigned, N> mk_inv(std::array<Index, N> const& neighs) {
  std::array<unsigned, N> ans;
  for (int i = 0; i < N; ++i) {
    for (int j = i; j < N; ++j) {
      if (neighs[i] == -neighs[j]) {
        ans[i] = j;
        ans[j] = i;
        break;
      }
    }
  }
  return ans;
}
}  // namespace detail

struct Neighbours {
  // shortcut to geometry class
  using gmy = hemelb::io::formats::geometry;
  using DoubleDisplacementArray =
      std::array<Vector, gmy::NumberOfDisplacements>;

  // The lattice vectors to neighbouring points (ordered as per the GMY file).
  static constexpr gmy::DisplacementArray Displacements = gmy::Neighbourhood;
  // Double version of the above
  static constexpr DoubleDisplacementArray DoubleDisplacements =
      detail::i2d(gmy::Neighbourhood);
  // Index of inverses, i.e. Displacement[i] == -Displacement[Inverses[i]]
  static constexpr std::array<unsigned, gmy::NumberOfDisplacements> Inverses =
      detail::mk_inv(gmy::Neighbourhood);
};

}  // namespace hemelb::gmytool::oct
#endif
