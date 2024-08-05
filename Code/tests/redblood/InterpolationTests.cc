// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "redblood/Interpolation.h"
#include "redblood/VelocityInterpolation.h"
#include "lb/lattices/D3Q15.h"
#include "lb/kernels/LBGK.h"
#include "lb/kernels/GuoForcingLBGK.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"

namespace hemelb::tests
{
    using namespace redblood;
    struct PlanarFunction
    {
        LatticePosition operator()(LatticePosition const &pos) const
        {
            return {
                    Dot(pos, util::Vector3D{1, 1, 1}),
                    Dot(pos, util::Vector3D{-1, 2, 1}),
                    Dot(pos, util::Vector3D{0, 0, 1})
            };
        }
        // We need to accept grid point and continuous input points
        LatticePosition operator()(LatticeVector const &pos) const
        {
            return operator()(pos.as<double>());
        }
        LatticePosition operator()(Dimensionless x, Dimensionless y, Dimensionless z) const
        {
            return operator()(LatticePosition(x, y, z));
        }
    };

    struct QuadraticFunction
    {
        // We need to accept grid point and continuous input points
        LatticePosition operator()(LatticePosition const &pos) const
        {
            Dimensionless const offset(0);
            auto dot_shift_square = [&](LatticePosition v) {
                auto tmp = Dot(pos, v) - offset;
                return tmp * tmp;
            };
            return {
                dot_shift_square({1, 1, 1}),
                dot_shift_square({0, 1, 0}),
                dot_shift_square({0, 0, 1})
            };
      }
      // We need to accept grid point and continuous input points
      LatticePosition operator()(LatticeVector const &pos) const
      {
	return operator()(pos.as<double>());
      }
      LatticePosition operator()(Dimensionless x, Dimensionless y, Dimensionless z) const
      {
	return operator()(LatticePosition(x, y, z));
      }
    };
    std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::FourPoint const &)
    {
      return {
              { LatticeVector(55, 51, 14), 0},
              { LatticeVector(55, 51, 15), 1},
              { LatticeVector(55, 51, 16), 1},
              { LatticeVector(55, 51, 17), 1},
              { LatticeVector(55, 52, 14), 1},
              { LatticeVector(56, 51, 14), 12},
              { LatticeVector(58, 54, 17), 47}
            };
    }
    std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::CosineApprox const &)
    {
      return offLatticeData(stencil::FourPoint());
    }
    std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::ThreePoint const &)
    {
      return
	{
              { LatticeVector(56, 52, 14), 0},
              { LatticeVector(56, 52, 15), 1},
              { LatticeVector(56, 52, 16), 1},
              { LatticeVector(56, 53, 14), 1},
              { LatticeVector(57, 52, 14), 6},
              { LatticeVector(58, 54, 16), 17}
	};
    }
    std::vector<std::pair<LatticeVector, size_t>> offLatticeData(stencil::TwoPoint const &)
    {
      return
	{
              { LatticeVector(56, 52, 15), 0},
              { LatticeVector(56, 52, 16), 1},
              { LatticeVector(56, 53, 15), 1},
              { LatticeVector(57, 52, 15), 2},
              { LatticeVector(57, 53, 16), 3}
	};
    }

    TEST_CASE("InterpolationTests", "[redblood]") {

      SECTION("testIndexIterator") {
	int constexpr n = 8;
	LatticeVector vectors[] = { LatticeVector(4, 3, 2),
				    LatticeVector(4, 3, 3),
				    LatticeVector(4, 3, 4),
				    LatticeVector(4, 4, 2),
				    LatticeVector(4, 5, 2),
				    LatticeVector(5, 3, 2),
				    LatticeVector(6, 3, 2),
				    LatticeVector(6, 5, 4), };
	size_t incs[] = { 0, 1, 1, 1, 3, 3, 9, 8};
	static_assert(sizeof(incs) == n*sizeof(size_t), "array size mismatch");
	// Checks iteration goes through correct sequence
	IndexIterator iterator(LatticeVector(5, 4, 3), 1);

	for (size_t i = 0; i < n; ++i) {
	  for (size_t j(0); j < incs[i]; ++j, ++iterator) {}

	  REQUIRE(LatticeVector::Zero() == (*iterator - vectors[i]));
	  REQUIRE(iterator.IsValid());
	}

	// Checks iterator becomes invalid
	++iterator;
	REQUIRE(not iterator.IsValid());
      }
    }

    auto approx = Approx(0).margin(1e-8);
    using StencilTypes = std::tuple<stencil::FourPoint, stencil::CosineApprox, stencil::ThreePoint, stencil::TwoPoint>;

    TEMPLATE_LIST_TEST_CASE("test Off Lattice interpolation", "[redblood]", StencilTypes) {
      using STENCIL = TestType;
      LatticePosition const pos(56.51, 52.9, 15.2);
      InterpolationIterator<STENCIL> iterator(pos);

      // Expected point and number of times to increment iterator
      auto const expected = offLatticeData(STENCIL());

      // Checks iteration goes through correct sequence
      for (auto const &item : expected)
	{
	  for (size_t j(0); j < item.second; ++j, ++iterator)
	    ;

	  REQUIRE(iterator.IsValid());
	  REQUIRE(item.first.x() == iterator->x());
	  REQUIRE(item.first.y() == iterator->y());
	  REQUIRE(item.first.z() == iterator->z());
	  REQUIRE(iterator.weight() == approx(STENCIL::stencil(pos - item.first)));
	}

      ++iterator;
      REQUIRE(not iterator.IsValid());
    }

    std::vector<LatticeVector> outsidePoints(stencil::FourPoint const&)
    {
      return {
	LatticeVector(57, 53, 13), LatticeVector(57, 53, 18), LatticeVector(57, 50, 17),
	LatticeVector(57, 55, 17), LatticeVector(54, 53, 17), LatticeVector(59, 53, 17)
      };
    }
    std::vector<LatticeVector> outsidePoints(stencil::CosineApprox const&)
    {
      return outsidePoints(stencil::FourPoint());
    }
    std::vector<LatticeVector> outsidePoints(stencil::ThreePoint const&)
    {
      return {
	LatticeVector(56, 52, 13), LatticeVector(56, 52, 17), LatticeVector(56, 51, 16),
	LatticeVector(56, 55, 16), LatticeVector(55, 52, 16), LatticeVector(59, 52, 16)
      };
    }
    std::vector<LatticeVector> outsidePoints(stencil::TwoPoint const&)
    {
      return {
	LatticeVector(56, 52, 14), LatticeVector(56, 52, 17), LatticeVector(56, 51, 15),
	LatticeVector(56, 54, 16), LatticeVector(55, 52, 16), LatticeVector(58, 52, 16)
      };
    }

    TEMPLATE_LIST_TEST_CASE("Interpolation tests OffLatticeZeroOutsideStencil", "[redblood]", StencilTypes)
    {
      using STENCIL = TestType;
      LatticePosition const pos(56.51, 52.9, 15.2);
      InterpolationIterator<STENCIL> iterator(pos);

      // Checks that outside iteration box, weights are zero
      for (auto const& vec : outsidePoints(STENCIL()))
	{
	  LatticeVector const dx(1, 0, 0), dy(0, 1, 0), dz(0, 0, 1);
	  REQUIRE(approx(0e0) == STENCIL::stencil(pos - vec));
	  // checks we are one step outside the iteration box only.
	  // this is really a test on the zero_vecs data, eg a test of the test.
	  size_t const one_non_zero = 0
	    + size_t(STENCIL::stencil(pos + dx - vec) != Approx(0))
	    + size_t(STENCIL::stencil(pos - dx - vec) != Approx(0))
	    + size_t(STENCIL::stencil(pos + dy - vec) != Approx(0))
	    + size_t(STENCIL::stencil(pos - dy - vec) != Approx(0))
	    + size_t(STENCIL::stencil(pos + dz - vec) != Approx(0))
	    + size_t(STENCIL::stencil(pos - dz - vec) != Approx(0));
	  REQUIRE(size_t(1) == one_non_zero);
	}
    }

    template<class FUNCTION, class STENCIL>
    void check(Dimensionless x, Dimensionless y, Dimensionless z,
	       Dimensionless tolerance = 1e-8)
    {
      FUNCTION func;
      LatticePosition expected(func(x, y, z));
      LatticePosition actual(interpolate<FUNCTION, STENCIL>(func, x, y, z));
      REQUIRE(Approx(expected.x()).margin(tolerance) == actual.x());
      REQUIRE(Approx(expected.y()).margin(tolerance) == actual.y());
      REQUIRE(Approx(expected.z()).margin(tolerance) == actual.z());
    }

    // Test interpolation when the point is on the grid
    TEMPLATE_LIST_TEST_CASE("InterpolateLinearFunction", "[redblood]", StencilTypes) {
      using STENCIL = TestType;

      SECTION("Linear function") {
	auto const tolerance = std::is_same<STENCIL, stencil::CosineApprox>::value ?
	  5e-2 :
	  1e-8;
	check<PlanarFunction, STENCIL>(0, 0, 0, tolerance);
	check<PlanarFunction, STENCIL>(0.1, 0.5, 0.6, tolerance);
	check<PlanarFunction, STENCIL>(-5.1, 0.5, 8.7, tolerance);
	check<PlanarFunction, STENCIL>(-5, 0, -1, tolerance);
      }
    
      SECTION("QuadraticFunction") {
	QuadraticFunction quad;
	// Error depends on variation on scale larger than stencil
	Dimensionless const tolerance( (quad(0, 0, 0) - quad(12, 12, 12)).GetMagnitude()
				       * 1e-2);
	check<QuadraticFunction, STENCIL>(0, 0, 0, tolerance);
	check<QuadraticFunction, STENCIL>(0.1, 0.5, 0.6, tolerance);
	check<QuadraticFunction, STENCIL>(-5.1, 0.5, 8.7, tolerance);
	check<QuadraticFunction, STENCIL>(-5, 0, -1, tolerance);
      }

      SECTION("MinMaxPosCornerCase") {
	// LatticeDistance 0x402effffffffffff is one bit short of 15.5. The
	// previous implementations of minimumPosImpl/maximumPosImpl were
	// returning a range of 4 lattice sites around it for the 3-point stencil.
	uint64_t const ux = 0x402effffffffffff;
	assert(sizeof(LatticeDistance) == sizeof(uint64_t));
	LatticeDistance const x(*reinterpret_cast<LatticeDistance const *>(&ux));

	auto const stencil_range = STENCIL::GetRange();
	auto const x_min = minimumPosImpl(x, stencil_range);
	auto const x_max = maximumPosImpl(x, stencil_range);
	decltype(stencil_range) const num_positions = x_max - x_min + 1;
	REQUIRE(stencil_range == num_positions);
      }
    }

    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<>, "VelocityInterpolationTests", "[redblood]") {
      
      using D3Q15 = lb::D3Q15;
      using LBGK = lb::LBGK<D3Q15>;
      using GuoForcingLBGK = lb::GuoForcingLBGK<D3Q15>;
      //     CPPUNIT_TEST_SUITE (VelocityInterpolationTests);
      //     CPPUNIT_TEST (testVelocityDataFromLatticeWithForces);
      //     CPPUNIT_TEST (testVelocityDataFromLatticeWithoutForces);CPPUNIT_TEST_SUITE_END();

      //   public:
      //     void setUp()
      //     {
      //       FourCubeBasedTestFixture::setUp();

      LatticeVector const min(dom->GetGlobalSiteMins());
      LatticeVector const max(dom->GetGlobalSiteMaxes());

      for (int dists = 0; dists < 2; ++dists) {
	for (site_t i = min[0]; i <= max[0]; ++i)
	  for (site_t j = min[1]; j <= max[1]; ++j)
	    for (site_t k = min[2]; k <= max[2]; ++k) {
	      distribn_t fold[D3Q15::NUMVECTORS];
	      std::fill(fold, fold+D3Q15::NUMVECTORS, 0.0);

	      fold[2] = i;
	      fold[4] = j;
	      fold[6] = k;
	      site_t const local = dom->GetContiguousSiteId(LatticeVector(i, j, k));
	      latDat->SetFOld<D3Q15>(local, fold);
	    }
	latDat->SwapOldAndNew();
      }
      
      //template<class KERNEL>
      auto velocityFromLatticeDataTester = [this](auto kernel, bool doforce)
      {
	using hemelb::redblood::details::VelocityFromLatticeData;
	using hemelb::redblood::details::HasForce;

	using KERNEL = decltype(kernel);
	// Check that type traits works as expected
	REQUIRE(HasForce<KERNEL>::value == doforce);

	//KERNEL kernel(initParams);
	VelocityFromLatticeData<KERNEL> velocityFunctor(*latDat);
    auto [begin, end] = dom->GetMidDomainSiteRange(0);
    auto const N = end - begin;

	REQUIRE(N > 0);

	for (auto index = begin; index < end; ++index) {
	  geometry::Site<geometry::FieldData const> const site(index, *latDat);

	  // Value to test
	  LatticeVelocity const actual = velocityFunctor(index);

	  // Check against directly computed values:
	  // We know how the lattice was setup
	  LatticeVector const position(site.GetGlobalSiteCoords());
	  LatticeVelocity const momentum = LatticeVelocity(D3Q15::CX[2],
							   D3Q15::CY[2],
							   D3Q15::CZ[2])
	    * Dimensionless(position[0])
	    + LatticeVelocity(D3Q15::CX[4], D3Q15::CY[4], D3Q15::CZ[4])
	    * Dimensionless(position[1])
	    + LatticeVelocity(D3Q15::CX[6], D3Q15::CY[6], D3Q15::CZ[6])
	    * Dimensionless(position[2]) + (doforce ?
					    site.GetForce() * 0.5 :
					    LatticeVelocity(0, 0, 0));
	  LatticeDensity const density(position[0] + position[1] + position[2]);

	  REQUIRE(actual == ApproxV(momentum / density));

	  // Check against Kernel + HydroVars implementation
	  typename KERNEL::VarsType hydroVars(site);
	  kernel.CalculateDensityMomentumFeq(hydroVars, site.GetIndex());

	  REQUIRE((actual ==  ApproxV(hydroVars.velocity)));
	}
      };

      SECTION("VelocityDataFromLatticeWithoutForces") {
	helpers::LatticeDataAccess(latDat.get()).ZeroOutForces();
	LBGK kernel{initParams};
	velocityFromLatticeDataTester(kernel, false);
      }
      SECTION("VelocityDataFromLatticeWithForces") {
            auto [begin, end] = dom->GetMidDomainSiteRange(0);

	for (size_t i = begin; i < end; ++i)
	  latDat->GetSite(i).SetForce(LatticeForceVector(i, 2 * i, double(i * i) * 0.0001));

	GuoForcingLBGK kernel{initParams};
	velocityFromLatticeDataTester(kernel, true);
      }
    }

}

