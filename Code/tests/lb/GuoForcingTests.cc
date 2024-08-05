// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <random>
#include <numeric>
#include <catch2/catch.hpp>

#include "lb/lattices/D3Q19.h"
#include "lb/lattices/Lattice.h"
#include "lb/Streamers.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "tests/helpers/FourCubeBasedTestFixture.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb::tests
{
    using fourcube =  helpers::FourCubeBasedTestFixture<>;
    class GuoForcingTests : public fourcube
    {
        using LatticeType = lb::D3Q19;
        using Kernel = lb::GuoForcingLBGK<LatticeType>;
        using VarsType = Kernel::VarsType;
        using Collision = lb::Normal<Kernel>;
        //using Lattice = lb::Lattice<LatticeType>;

        std::unique_ptr<lb::MacroscopicPropertyCache> propertyCache;
        Approx apprx = Approx(0.0).margin(1e-8);
        ApproxVector<double> apprx_vec = ApproxVector<double>{0.0}.Margin(1e-8);

    public:
        GuoForcingTests() : fourcube{}
        {
            // Set forces equal to site index
            size_t const nFluidSites = dom->GetLocalFluidSiteCount();
            for (size_t site(0); site < nFluidSites; ++site) {
                latDat->GetSite(site).SetForce(GetMeAForce(site));
            }
            propertyCache = std::make_unique<lb::MacroscopicPropertyCache>(*simState, *dom);
        }

        // Tests the correct force is assigned to HydroVars, and that the right
        // specialization is used as well.
        void testHydroVarsGetForce()
        {
            size_t const nFluidSites = dom->GetLocalFluidSiteCount();
            for (size_t i(1); i < nFluidSites; i <<= 1) {
                // REQUIRE(helpers::is_zero(HydroVars(latDat->GetSite(i)).force - GetMeAForce(i)));
                REQUIRE(VarsType(latDat->GetSite(i)).force == apprx_vec(GetMeAForce(i)));
            }
        }

        // Checks equation 19 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
        void testCalculatDensityAndMomentum()
        {
            // Create fake density and compute resulting momentum and velocity
            // velocity is not yet corrected for forces
            auto const [f, expected_density, expected_momentum] = []() {
                typename LatticeType::FArray f;
                PhysicalVelocity expected_momentum(0, 0, 0);
                PhysicalDensity expected_density(0);
                for (size_t i(0); i < LatticeType::NUMVECTORS; ++i) {
                    f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
                    expected_density += f[i];
                    expected_momentum += LatticeType::VECTORS[i].as<PhysicalSpeed>() * f[i];
                }
                return std::make_tuple(f, expected_density, expected_momentum);
            }();

            PhysicalVelocity const force(1, 2, 3);
            PhysicalVelocity momentum(0, 0, 0);
            PhysicalDensity density(0);

            // Check when forces are zero
            LatticeType::CalculateDensityAndMomentum(f,
                                                     LatticeForceVector::Zero(),
                                                     density,
                                                     momentum);
            REQUIRE(density == apprx(expected_density));
            REQUIRE(momentum == apprx_vec(expected_momentum));

            // Check when forces are not zero
            LatticeType::CalculateDensityAndMomentum(f,
                                                     force,
                                                     density,
                                                     momentum);
            REQUIRE(density == apprx(expected_density));
            REQUIRE(momentum == apprx_vec(expected_momentum + force * 0.5));
        }

        // Checks Guo Correction is included in FEq and velocity
        void testCalculateDensityMomentumFEq()
        {
            LatticeForceVector const force(1, 2, 3);
            auto const f = [&]() {
                typename LatticeType::FArray f;
                for (size_t i = 0; i < LatticeType::NUMVECTORS; ++i)
                    f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
                return f;
            }();

            PhysicalVelocity noforce_momentum(0, 0, 0), force_momentum(0, 0, 0);
            PhysicalVelocity noforce_velocity(0, 0, 0), force_velocity(0, 0, 0);
            PhysicalDensity noforce_density(0), force_density(0);
            distribn_t noforce_feq[LatticeType::NUMVECTORS], force_feq[LatticeType::NUMVECTORS],
                    feq[LatticeType::NUMVECTORS];

            // Compute without force
            LatticeType::CalculateDensityMomentumFEq(f,
                                                     noforce_density,
                                                     noforce_momentum,
                                                     noforce_velocity,
                                                     noforce_feq);
            // Compute with forces
            LatticeType::CalculateDensityMomentumFEq(f,
                                                     force,
                                                     force_density,
                                                     force_momentum,
                                                     force_velocity,
                                                     force_feq);
            REQUIRE(noforce_density == apprx(force_density));
            REQUIRE(noforce_momentum + force * 0.5 == apprx_vec(force_momentum).Margin(1e-8));
            REQUIRE(apprx_vec(force_velocity) == force_momentum / force_density);

            // Compute feq with forces directly
            LatticeType::CalculateFeq(force_density,
                                      force_momentum,
                                      feq);
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
                REQUIRE(force_feq[i] == apprx(feq[i]));
        }

        void ForceDistributionTestInstance(LatticeVelocity const &_velocity, const distribn_t _tau,
                                           const LatticeForceVector &_force)
        {
            distribn_t expected_Fi[LatticeType::NUMVECTORS];
            distribn_t actual_Fi[LatticeType::NUMVECTORS];

            ForceDistribution(_velocity, _tau, _force, expected_Fi);
            LatticeType::CalculateForceDistribution(_tau,
                                                    _velocity,
                                                    _force,
                                                    actual_Fi);

            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
                REQUIRE(actual_Fi[i] == apprx(expected_Fi[i]));
        }

        // Equation 20 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
        void testForceDistributionRegression()
        {
            ForceDistributionTestInstance(LatticeVelocity(1, 0, 0),
                                          0.25,
                                          LatticeForceVector(1.5, 0, 0));
            ForceDistributionTestInstance(LatticeVelocity(0, 1, 0),
                                          0.25,
                                          LatticeForceVector(1.5, 0, 0));
            ForceDistributionTestInstance(LatticeVelocity(0, 0, 1),
                                          0.25,
                                          LatticeForceVector(3.0, 0, 0));
        }

        // tau = 0.5 => zero force: look at prefactor
        void testForceDistributionBadTau()
        {
            distribn_t Fi[LatticeType::NUMVECTORS];
            LatticeType::CalculateForceDistribution(0.5, {1.0, 10.0, 100.0}, {1, 1, 1}, Fi);
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
                REQUIRE(Fi[i] == apprx(0));
        }

        // velocity and force colinear with CX, CY, CZ:
        // tests second term in eq 20
        void testForceDistributionColinear()
        {
            distribn_t Fi[LatticeType::NUMVECTORS];
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i) {
                LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
                LatticeType::CalculateForceDistribution(0.25, ei, ei, Fi);
                distribn_t ei_norm(Dot(ei, ei));
                distribn_t const expected = (1. - 0.5 / 0.25) * LatticeType::EQMWEIGHTS[i]
                                            * (ei_norm * ei_norm * 9.0);
                REQUIRE(Fi[i] == apprx(expected));
            }
        }
        // force colinear with CX, CY, CZ, and zero force:
        // tests first term in eq 20
        void testForceDistributionZeroVelocity()
        {
            distribn_t Fi[LatticeType::NUMVECTORS];
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i) {
                LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
                LatticeType::CalculateForceDistribution(0.25, LatticeVelocity::Zero(), ei, Fi);
                distribn_t ei_norm(Dot(ei, ei));
                distribn_t const expected = (1. - 0.5 / 0.25) * LatticeType::EQMWEIGHTS[i] * (ei_norm * 3.0);
                REQUIRE(Fi[i] == apprx(expected));
            }
        }

        void testNoForceMeansNoForceDistribution()
        {
            auto const N = 10;
            std::random_device rd;
            std::default_random_engine e1(rd());
            std::uniform_real_distribution<distribn_t> rdm(-10, 10);

            distribn_t Fi[LatticeType::NUMVECTORS];
            for(size_t i = 0; i < N; ++i) {
	      LatticeType::CalculateForceDistribution(1e0, {rdm(e1), rdm(e1), rdm(e1)}, LatticeForceVector::Zero(), Fi);
                for(size_t j(0); j < LatticeType::NUMVECTORS; ++j) {
                    REQUIRE(Fi[j] == Approx(0).margin(1e-12));
                }
            }
        }

        void testZerothOrderMomentOfForceIsNull()
        {
            auto const N = 10;
            std::random_device rd;
            std::default_random_engine e1(rd());
            std::uniform_real_distribution<distribn_t> rdm(-10, 10);

            distribn_t Fi[LatticeType::NUMVECTORS];
            // No velocity
            for (size_t i = 0; i < N; ++i) {
                LatticeType::CalculateForceDistribution(
		    1e0, LatticeVelocity::Zero(),
		    {rdm(e1), rdm(e1), rdm(e1)},
		    Fi);
                auto const moment = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
                REQUIRE(moment == apprx(0.0));
            }
            // Force and Velocity are perpendicular
            for(size_t i = 0; i < N; ++i) {
                auto const a = rdm(e1), b = rdm(e1);
                LatticeType::CalculateForceDistribution(1e0, {rdm(e1), a, b}, {0, -b, a}, Fi);
                auto const moment0 = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
                REQUIRE(moment0 == apprx(0.0));

                LatticeType::CalculateForceDistribution(1e0, {0, a, b}, {rdm(e1), -b, a}, Fi);
                auto const moment1 = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
                REQUIRE(moment1 == apprx(0.0));
            }

            // Force and velocity tied via pops
            // populations
            for (size_t i = 0; i < 1; ++i) {
                LatticeForceVector const F(rdm(e1), rdm(e1), rdm(e1));
                LatticeVelocity v = F * 0.5;
                LatticeDensity rho(0);
                for(size_t j(0); j < LatticeType::NUMVECTORS; ++j) {
                    auto const fi = rdm(e1);
                    rho += fi;
                    v += LatticeVelocity(LatticeType::CX[j], LatticeType::CY[j], LatticeType::CZ[j]) * fi;
                }
                v = v / rho;
                LatticeType::CalculateForceDistribution(1e0, v, F, Fi);
                auto const moment = std::accumulate(Fi, Fi+LatticeType::NUMVECTORS, 0e0);
                REQUIRE(moment == apprx(0.0));
            }
        }

        void testGuoEquation()
        {
            std::random_device rd;
            std::default_random_engine e1(rd());
            std::uniform_real_distribution<distribn_t> rdm(-10, 10);

            LatticeVelocity const v(rdm(e1), rdm(e1), rdm(e1));
            LatticeForceVector const F(rdm(e1), rdm(e1), rdm(e1));
            LatticeForce sumFi(0);
            for(size_t j = 0; j < LatticeType::NUMVECTORS; ++j) {
                LatticeVelocity const ei(LatticeType::CXD[j], LatticeType::CYD[j], LatticeType::CZD[j]);
                sumFi += LatticeType::EQMWEIGHTS[j] * Dot((ei - v) * 3e0 + ei * Dot(ei, v) * 9e0, F);
            }
            REQUIRE(sumFi == apprx(0.0));
        }

        // Checks collision process includes forcing
        // Assumes that equilibrium and non-equilibrium populations, velocity,
        // and forces are computed elsewhere.
        // This just makes sure the force distribution included correctly.
        void testDoCollide()
        {
            distribn_t f[LatticeType::NUMVECTORS], Fi[LatticeType::NUMVECTORS];
            LatticeForceVector const force(1, 4, 8);
            VarsType hydroVars(f, force);
            for (size_t i = 0; i < LatticeType::NUMVECTORS; ++i) {
                f[i] = distribn_t(std::rand()) / distribn_t(RAND_MAX);
                hydroVars.SetFNeq(i, distribn_t(std::rand()) / distribn_t(RAND_MAX));
            }

            hydroVars.velocity = LatticeVelocity(1, 2, 3);
            ForceDistribution(hydroVars.velocity, lbmParams.GetTau(), hydroVars.force, Fi);

            Kernel kernel(initParams);
            kernel.Collide(&lbmParams, hydroVars);
            for (size_t i = 0; i < LatticeType::NUMVECTORS; ++i) {
                REQUIRE(
                        hydroVars.GetFPostCollision()[i] ==
                        apprx(f[i] +
                              hydroVars.GetFNeq()[i] * lbmParams.GetOmega() +
                              Fi[i])
                );
            }
        }

        void testSimpleCollideAndStream()
        {
            // Lattice with a single non-zero distribution in the middle
            // Same for forces
            LatticeVector const position(2, 2, 2);
            auto&& site = latDat->GetSite(position);
            auto i_site = site.GetIndex();
            helpers::allZeroButOne<LatticeType>(*latDat, position);

            // Get collided distributions at this site
            // Will compare zero and non-zero forces, to make sure they are different
            // Assumes collides works, since tested in TestDoCollide
            LatticeType::FArray withForce, withoutForce;
            FPostCollision(site.GetFOld<LatticeType>(), site.GetForce(), withForce);
            FPostCollision(site.GetFOld<LatticeType>(), LatticeForceVector(0, 0, 0), withoutForce);

            // Stream that site
            using lb::BulkStreamer;
            BulkStreamer<Collision> streamer(initParams);
            streamer.StreamAndCollide(i_site, i_site + 1, &lbmParams, *latDat, *propertyCache);

            // Now check streaming worked correctly
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i) {
                auto const index = site.GetStreamedIndex<LatticeType>(i);
                REQUIRE(withForce[i] == apprx(*helpers::GetFNew(*latDat, index)));
                // And that forces from streaming site were used
                REQUIRE(withForce[i] != apprx(withoutForce[i]));
            }
        }

        void testSimpleBounceBack()
        {
            // Lattice with a single non-zero distribution in the middle
            // Same for forces
            LatticeVector const position(1, 1, 1);
            auto&& site = latDat->GetSite(position);
            auto i_site = site.GetIndex();
            helpers::allZeroButOne<LatticeType>(*latDat, position);

            // Get collided distributions at this site
            // Will compare zero and non-zero forces, to make sure they are different
            // Assumes collides works, since tested in TestDoCollide
            distribn_t withForce[LatticeType::NUMVECTORS], withoutForce[LatticeType::NUMVECTORS];
            FPostCollision(site.GetFOld<LatticeType>(), site.GetForce(), withForce);
            FPostCollision(site.GetFOld<LatticeType>(), LatticeForceVector(0, 0, 0), withoutForce);

            // Stream that site
            using SBB = lb::StreamerTypeFactory<
                    lb::BounceBackLink<Collision>,
                    lb::NullLink<Collision>
            >;
            SBB streamer(initParams);
            streamer.StreamAndCollide(i_site, i_site + 1, &lbmParams, *latDat, *propertyCache);

            distribn_t const * const actual = helpers::GetFNew<LatticeType>(*latDat, position);
            bool paranoia = false;
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i) {
                if (not site.HasWall(i))
                    continue;
                paranoia = true;
                REQUIRE(withForce[i] == apprx(actual[LatticeType::INVERSEDIRECTIONS[i]]));
                // And that forces from streaming site were used
                REQUIRE(withForce[i] != apprx(withoutForce[i]));
            }
            REQUIRE(paranoia);
        }

    protected:
        LatticeForceVector GetMeAForce(size_t site)
        {
            return LatticeForceVector(LatticeForce(site),
                                      LatticeForce(site) + 0.1,
                                      LatticeForce(site) + 0.3);
        }

        // Sets up _force distribution from actual 3d vector
        void ForceDistribution(LatticeVelocity const &_velocity, const distribn_t _tau,
                               const LatticeForceVector &_force, distribn_t Fi[])
        {
            const distribn_t inv_cs2(3e0);
            const distribn_t inv_cs4(9e0);
            const distribn_t prefactor(1. - 0.5 / _tau);
            for (size_t i(0); i < LatticeType::NUMVECTORS; ++i)
            {
                LatticeVelocity const ei(LatticeType::CX[i], LatticeType::CY[i], LatticeType::CZ[i]);
                LatticeForceVector const forcing( (ei - _velocity) * inv_cs2
                                                  + ei * (Dot(ei, _velocity) * inv_cs4));
                Fi[i] = Dot(forcing, _force) * LatticeType::EQMWEIGHTS[i] * prefactor;
            }
        }

        void FPostCollision(LatticeType::const_span _f, LatticeForceVector const &_force,
                            LatticeType::mut_span _fout)
        {
            VarsType hydroVars(_f, _force);
            Kernel kernel(initParams);
            // Index is never used ... so give a fake value
            kernel.CalculateDensityMomentumFeq(hydroVars, 999999);
            kernel.Collide(&lbmParams, hydroVars);

            std::copy(hydroVars.GetFPostCollision().begin(),
                      hydroVars.GetFPostCollision().end(),
                      _fout.begin());
        }

    };

    METHOD_AS_TEST_CASE(GuoForcingTests::testHydroVarsGetForce,
			"Tests the correct force is assigned to HydroVars, and that the right specialization is used as well.",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testCalculatDensityAndMomentum,
			"Checks equation 19 of Guo paper (doi: 10.1103/PhysRevE.65.046308)",
			"[redblood]");

    METHOD_AS_TEST_CASE(GuoForcingTests::testCalculateDensityMomentumFEq,
			"Checks Guo Correction is included in FEq and velocity",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testForceDistributionRegression,
			"Equation 20 of Guo paper (doi: 10.1103/PhysRevE.65.046308)",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testForceDistributionBadTau,
			"tau = 0.5 => zero force: look at prefactor",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testForceDistributionColinear,
			"velocity and force colinear with CX, CY, CZ: tests second term in eq 20",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testForceDistributionZeroVelocity,
			"force colinear with CX, CY, CZ, and zero force: tests first term in eq 20",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testZerothOrderMomentOfForceIsNull,
			"test zeroth order moment of force is zero",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testDoCollide,
			"Checks collision process includes forcing",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testNoForceMeansNoForceDistribution,
			"test that no force means no contribution to force in distributions",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testGuoEquation,
			"test Guo equation",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testSimpleCollideAndStream,
			"test that the addition of streaming works",
			"[redblood]");
    METHOD_AS_TEST_CASE(GuoForcingTests::testSimpleBounceBack,
			"test that works in the presence of a wall + SBB",
			"[redblood]");
}
