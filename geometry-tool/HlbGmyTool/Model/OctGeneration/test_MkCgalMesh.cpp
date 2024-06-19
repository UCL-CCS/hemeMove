// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "MkCgalMesh.h"
#include "TestResources/Meshes.h"

namespace hemelb::gmytool::oct {

TEST_CASE("SphereCgalMesh") {
  auto sphere = SimpleMeshFactory::MkSphere();
  auto surface = MkCgalMesh(sphere->points, sphere->triangles);
  REQUIRE(surface->is_closed());
  REQUIRE(surface->is_pure_triangle());
  REQUIRE(surface->is_valid());
}

}  // namespace hemelb::gmytool::oct
