// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/RBCInserter.h"
#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb
{
  namespace tests
  {

    TEST_CASE_METHOD(helpers::HasCommsTestFixture,
		     "RBCInserterWithPerturbationParallelTests",
		     "[redblood]") {
      // This test ensures that all the processes in a simulation add the same
      // amount of delta_x and delta_y to any given cell before trying to drop it.
      using namespace hemelb::redblood;

      auto const deltaX = LatticePosition(2, 0, 0);
      auto const deltaY = LatticePosition(0, 4, 0);
      auto const identity = rotationMatrix(LatticePosition(0, 0, 1),
                                                 LatticePosition(0, 0, 1));

      TemplateCellContainer cells;
      cells.emplace("joe", std::make_shared<Cell>(tetrahedron()));
      RBCInserterWithPerturbation inserter([](){ return true;}, cells["joe"]->clone(),
					   identity, 0e0, 0e0, deltaX, deltaY);

      // Original implementation of drop() used rand(), which was generating different
      // numbers across processes in drop(), even when seeded equally at the beginning
      // of the simulation across processes. This was probably due to rand() being
      // called a different number of times before reaching drop() across processes.
      // Force pseudo-random number generator to be out of sync to test new random
      // number generation approach.
      auto const world = Comms();
      if (world.Rank() == 0)
	{
	  auto volatile r = rand();
	}

      auto const cell = inserter.drop();
      auto const barycentre = cell->GetBarycentre();

      const double p = world.Size();
      auto const x_times_p = world.AllReduce(barycentre.x(), MPI_SUM);
      auto const y_times_p = world.AllReduce(barycentre.y(), MPI_SUM);
      auto const z_times_p = world.AllReduce(barycentre.z(), MPI_SUM);
      auto const mean_barycentre = LatticePosition{x_times_p, y_times_p, z_times_p} / p;
      REQUIRE(barycentre == mean_barycentre);
    }
  }
}

