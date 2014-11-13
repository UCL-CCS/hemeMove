//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_GUOFORCING_H
#define HEMELB_UNITTESTS_GUOFORCING_H

#include <cppunit/TestFixture.h>
#include "lb/lattices/D3Q15.h"
#include "lb/lattices/Lattice.h"
#include "lb/kernels/GuoForcingLBGK.h"
#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "unittests/helpers/Comparisons.h"

namespace hemelb { namespace unittests {

  class GuoForcingTests
        : public helpers::FourCubeBasedTestFixture {
    CPPUNIT_TEST_SUITE(GuoForcingTests);
      CPPUNIT_TEST(testHydroVarsGetForce);
      CPPUNIT_TEST(testCalculatDensityAndMomentum);
      CPPUNIT_TEST(testCalculateDensityMomentumFEq);
    CPPUNIT_TEST_SUITE_END();
    public:
      void setUp() {
        FourCubeBasedTestFixture::setUp();
        // Set forces equal to site index
        size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
        for(size_t site(0); site < nFluidSites; ++site)
          latDat->GetSite(site).SetForce( GetMeAForce(site) );
      }

      // Tests the correct force is assigned to HydroVars, and that the right
      // specialization is used as well.
      void testHydroVarsGetForce() {
        size_t const nFluidSites(latDat->GetLocalFluidSiteCount());
        typedef lb::kernels::GuoForcingLBGK<lb::lattices::D3Q15> t_Kernel;
        typedef t_Kernel::KHydroVars t_HydroVars;
        for(size_t i(1); i < nFluidSites; i <<= 1)
          CPPUNIT_ASSERT(helpers::is_zero(
            t_HydroVars(latDat->GetSite(i)).force - GetMeAForce(i)
          ));
      }

      // Checks equation 19 of Guo paper (doi: 10.1103/PhysRevE.65.046308)
      void testCalculatDensityAndMomentum() {
        typedef lb::lattices::D3Q15 D3Q15;
        typedef lb::lattices::Lattice<D3Q15> Lattice;

        // Create fake density and compute resulting momentum and velocity
        // velocity is not yet corrected for forces
        distribn_t f[D3Q15::NUMVECTORS];
        PhysicalVelocity expected_momentum(0, 0, 0);
        PhysicalDensity expected_density(0);
        for(size_t i(0); i < lb::lattices::D3Q15::NUMVECTORS; ++i) {
          f[i] = rand();
          expected_density += f[i];
          expected_momentum +=
            PhysicalVelocity(D3Q15::CX[i], D3Q15::CY[i], D3Q15::CZ[i]) * f[i];
        }
        // ensures that functions take const value
        distribn_t const * const const_f = f;

        PhysicalVelocity const force(1, 2, 3);
        PhysicalVelocity momentum(0, 0, 0);
        PhysicalDensity density(0);

        // Check when forces are zero
        Lattice :: CalculateDensityAndMomentum(const_f, 0, 0, 0, density,
            momentum[0], momentum[1], momentum[2]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(density, expected_density, 1e-8);
        CPPUNIT_ASSERT(helpers::is_zero(momentum - expected_momentum));

        // Check when forces are not zero
        Lattice :: CalculateDensityAndMomentum(const_f, force[0], force[1],
            force[2], density, momentum[0], momentum[1], momentum[2]);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(density, expected_density, 1e-8);
        CPPUNIT_ASSERT(helpers::is_zero(
              momentum - expected_momentum - force * 0.5));
      }

      // Checks Guo Correction is included in FEq and velocity
      void testCalculateDensityMomentumFEq() {
        typedef lb::lattices::D3Q15 D3Q15;
        typedef lb::lattices::Lattice<D3Q15> Lattice;

        LatticeForceVector const force(1, 2, 3);
        distribn_t f[D3Q15::NUMVECTORS];
        for(size_t i(0); i < lb::lattices::D3Q15::NUMVECTORS; ++i)
          f[i] = rand();
        // ensures that functions take const value
        distribn_t const * const const_f = f;

        PhysicalVelocity noforce_momentum(0, 0, 0), force_momentum(0, 0, 0);
        PhysicalVelocity noforce_velocity(0, 0, 0), force_velocity(0, 0, 0);
        PhysicalDensity noforce_density(0), force_density(0);
        distribn_t noforce_feq[D3Q15::NUMVECTORS],
          force_feq[D3Q15::NUMVECTORS], feq[D3Q15::NUMVECTORS];

        // Compute without force
        Lattice :: CalculateDensityMomentumFEq(
            const_f,
            noforce_density,
            noforce_momentum[0], noforce_momentum[1], noforce_momentum[2],
            noforce_velocity[0], noforce_velocity[1], noforce_velocity[2],
            noforce_feq
        );
        // Compute with forces
        Lattice :: CalculateDensityMomentumFEq(
            const_f, force[0], force[1], force[2],
            force_density,
            force_momentum[0], force_momentum[1], force_momentum[2],
            force_velocity[0], force_velocity[1], force_velocity[2],
            force_feq
        );
        CPPUNIT_ASSERT_DOUBLES_EQUAL(noforce_density, force_density, 1e-8);
        CPPUNIT_ASSERT(helpers::is_zero(
          noforce_momentum + force * 0.5 - force_momentum
        ));
        CPPUNIT_ASSERT(helpers::is_zero(
          force_velocity - force_momentum / force_density
        ));

        // Compute feq with forces directly
        Lattice :: CalculateFeq(
            force_density,
            force_momentum[0], force_momentum[1], force_momentum[2],
            feq
        );
        for(size_t i(0); i < D3Q15::NUMVECTORS; ++i)
          CPPUNIT_ASSERT_DOUBLES_EQUAL(force_feq[i], feq[i], 1e-8);
      }

    protected:
      LatticeForceVector GetMeAForce(size_t site) {
        return LatticeForceVector(
            LatticeForce(site),
            LatticeForce(site) + 0.1,
            LatticeForce(site) + 0.3
        );
      }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION(GuoForcingTests);
}}

#endif
