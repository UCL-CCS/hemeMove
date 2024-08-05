// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "lb/IncompressibilityChecker.hpp"

#include "tests/lb/BroadcastMocks.h"
#include "tests/lb/LbTestsHelper.h"
#include "tests/reporting/Mocks.h"
#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb::tests
{
    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "IncompressibilityCheckerTests") {
      using LATTICE = lb::D3Q15;

      LbTestsHelper::InitialiseAnisotropicTestData<LATTICE>(*latDat);

      auto cache = std::make_unique<lb::MacroscopicPropertyCache>(*simState, *dom);
      cache->densityCache.SetRefreshFlag();
      cache->velocityCache.SetRefreshFlag();
      LbTestsHelper::UpdatePropertyCache<LATTICE>(*latDat, *cache, *simState);

      // These are the smallest and largest density values in
      // FourCubeLatticeData by default. The lattice class below must
      // be consistent with the one used in
      // FourCubeLatticeData.
      // 
      // TODO: Consider templating FourCubeLatticeData over lattice
      // class, so both can be controlled from the test.
      auto numDirections = distribn_t(LATTICE::NUMVECTORS);
      auto numSites = distribn_t(dom->GetLocalFluidSiteCount());

      // = sum_{j=1}^{numDirections} j/10 = 12 with current configuration of FourCubeLatticeData
      distribn_t smallestDefaultDensity = numDirections * (numDirections + 1) / 20;
      // = sum_{j=1}^{numDirections} j/10 + (numSites-1)/100 = 21.45 with current configuration of FourCubeLatticeData
      distribn_t largestDefaultDensity = (numDirections * (numDirections + 1) / 20) + ( (numSites - 1) * numDirections / 100);
      distribn_t largestDefaultVelocityMagnitude = 0.0433012701892219; // with current configuration of FourCubeLatticeData
      auto apprx = [](double x) {
	return Approx(x).margin(1e-9);
      };

      auto timings = std::make_unique<hemelb::reporting::Timers>();
      auto net = std::make_unique<net::Net>(Comms());

      auto AdvanceActorOneTimeStep = [&](net::IteratedAction& actor) {
	cache->densityCache.SetRefreshFlag();
	LbTestsHelper::UpdatePropertyCache<LATTICE>(*latDat, *cache, *simState);

	actor.RequestComms();
	actor.PreSend();
	actor.PreReceive();
	actor.PostReceive();
	actor.EndIteration();
      };
      
      SECTION("IncompressibilityCheckerRootNode") {
	lb::IncompressibilityChecker<net::BroadcastMockRootNode> incompChecker(dom,
									       net.get(),
									       simState.get(),
									       *cache,
									       *timings,
									       10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

	// Not available until the first broadcast has finished
	REQUIRE(!incompChecker.AreDensitiesAvailable());
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(smallestDefaultDensity) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(largestDefaultDensity) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx( (largestDefaultDensity - smallestDefaultDensity)
		       / hemelb::lb::REFERENCE_DENSITY) ==
		incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(largestDefaultVelocityMagnitude) == incompChecker.GetGlobalLargestVelocityMagnitude());

	// The broadcast mock injects some smaller and larger densities (1,100) coming from one of the children
	REQUIRE(incompChecker.AreDensitiesAvailable());
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(1.0) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(100.0) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx(99.0) == incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(!incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(10.0) == incompChecker.GetGlobalLargestVelocityMagnitude());

	// The previous values are not reported by any children anymore. Testing that the checker remembers them
	REQUIRE(incompChecker.AreDensitiesAvailable());
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(1.0) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(100.0) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx(99.0) == incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(!incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(10.0) == incompChecker.GetGlobalLargestVelocityMagnitude());
      }

      SECTION("IncompressibilityCheckerLeafNode") {
	lb::IncompressibilityChecker<net::BroadcastMockLeafNode> incompChecker(dom,
									       net.get(),
									       simState.get(),
									       *cache,
									       *timings,
									       10.0); // Will accept a max/min of (21.45, 12) but not (100,1)

	// First pass down the tree (uninitialised values being broadcasted)
	AdvanceActorOneTimeStep(incompChecker);

	// First pass up the tree (root node gets min/max values reported by the leaf node being simulated)
	AdvanceActorOneTimeStep(incompChecker);

	// Second pass down the tree (all nodes get the min/max values reported by the leaf node being simulated)
	AdvanceActorOneTimeStep(incompChecker);

	// This time the current leaf node reported global minimum and maximum values
	REQUIRE(apprx(smallestDefaultDensity) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(largestDefaultDensity) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx( (largestDefaultDensity - smallestDefaultDensity)
		       / hemelb::lb::REFERENCE_DENSITY) ==
				      incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(largestDefaultVelocityMagnitude) ==
				     incompChecker.GetGlobalLargestVelocityMagnitude());

	// Second pass up. The broadcast mock injects some smaller and larger densities (1,100) coming from another hypothetical leaf node.
	AdvanceActorOneTimeStep(incompChecker);

	// Third pass down. (1,100) arrives to all the leaf nodes.
	AdvanceActorOneTimeStep(incompChecker);

	REQUIRE(apprx(1.0) == incompChecker.GetGlobalSmallestDensity());
	REQUIRE(apprx(100.0) == incompChecker.GetGlobalLargestDensity());
	REQUIRE(apprx(99.0) == incompChecker.GetMaxRelativeDensityDifference());
	REQUIRE(!incompChecker.IsDensityDiffWithinRange());
	REQUIRE(apprx(10.0) == incompChecker.GetGlobalLargestVelocityMagnitude());
      }

    }
}
