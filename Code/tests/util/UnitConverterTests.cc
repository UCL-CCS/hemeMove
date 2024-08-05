// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "tests/lb/LbTestsHelper.h"

#include "util/UnitConverter.h"
#include "util/Matrix3D.h"
#include "constants.h"
#include "lb/lattices/Lattice.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb::tests
{
    using namespace hemelb::util;

    TEST_CASE("UnitConverterTests") {
      auto unitConverter = util::UnitConverter(1., 1., Vector3D<double>(0.), 1000.0, 0.0);
      PhysicalPressure pressPa = 81.0;
      LatticeDensity densityLatt = unitConverter.ConvertPressureToLatticeUnits(pressPa) / Cs2;
      distribn_t tau = 0.5;
      const double epsilon = 1e-9;
      auto apprx = [&](double x) {
	return Approx(x).margin(epsilon);
      };

      SECTION("TestPressure") {
	REQUIRE(apprx(pressPa) ==
		unitConverter.ConvertPressureToPhysicalUnits(densityLatt * Cs2));
      }

      SECTION("TestSimpleStressTensor") {
	auto fNonEquilibrium = LbTestsHelper::ZeroArray<lb::D3Q15>();

	util::Matrix3D stressTensor;
	lb::D3Q15::CalculateStressTensor(densityLatt, tau, fNonEquilibrium, stressTensor);
	util::Matrix3D stressTensorPhys = unitConverter.ConvertFullStressTensorToPhysicalUnits(stressTensor);

	REQUIRE(apprx(pressPa) == stressTensorPhys[0][0]);
	REQUIRE(apprx(pressPa) == stressTensorPhys[1][1]);
	REQUIRE(apprx(pressPa) == stressTensorPhys[2][2]);
	REQUIRE(apprx((distribn_t) 0.) == stressTensorPhys[1][0]);
      }

      SECTION("TestSimpleTractionVector") {
    auto fNonEquilibrium = LbTestsHelper::ZeroArray<lb::D3Q15>();
	util::Vector3D<Dimensionless> wallNormal(0.0);
	wallNormal[0] = 1.0;

	util::Vector3D<LatticeStress> traction;
	lb::D3Q15::CalculateTractionOnAPoint(densityLatt,
						       tau,
						       fNonEquilibrium,
						       wallNormal,
						       traction);
	util::Vector3D<PhysicalStress> tractionPhys = unitConverter.ConvertTractionToPhysicalUnits(traction,
												   wallNormal);

	REQUIRE(apprx(pressPa) == tractionPhys[0]);
	REQUIRE(apprx(0.0) == tractionPhys[1]);
	REQUIRE(apprx(0.0) == tractionPhys[2]);
      }

    }
}
