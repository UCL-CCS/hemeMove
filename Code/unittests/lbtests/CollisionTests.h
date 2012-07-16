// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_UNITTESTS_LBTESTS_COLLISIONTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_COLLISIONTESTS_H

#include <cppunit/TestFixture.h>


#include "unittests/helpers/FourCubeBasedTestFixture.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      /**
       * Class to test the collision operators. These tests are for the functions involved in
       * calculating the post-collision values, specifically CalculatePreCollision and Collide.
       * For each collision operator, we test that these functions return the expected values of
       * hydrodynamic quantities and f-distribution values.
       *
       * Note that we are only testing collision operators here, so we
       * can assume that the kernel objects work perfectly.
       */
      class CollisionTests : public helpers::FourCubeBasedTestFixture
      {
          CPPUNIT_TEST_SUITE( CollisionTests);
          CPPUNIT_TEST( TestNonZeroVelocityEquilibriumFixedDensity);
          CPPUNIT_TEST( TestZeroVelocityEquilibriumFixedDensity);
          CPPUNIT_TEST( TestZeroVelocityEquilibrium);
          CPPUNIT_TEST( TestNormal);CPPUNIT_TEST_SUITE_END();
        public:
          void setUp()
          {

            FourCubeBasedTestFixture::setUp();
            // Create the inlet and outlet boundary objects.
            inletBoundary = new lb::boundaries::BoundaryValues(geometry::INLET_TYPE,
                                                               latDat,
                                                               simConfig->GetInlets(),
                                                               simState,
                                                               unitConverter);
            outletBoundary = new lb::boundaries::BoundaryValues(geometry::OUTLET_TYPE,
                                                                latDat,
                                                                simConfig->GetOutlets(),
                                                                simState,
                                                                unitConverter);



            lbgk = new lb::kernels::LBGK<lb::lattices::D3Q15>(initParams);

            // Initialise all 4 types of conditions, using boundary objects for the collision types
            // that will need them.
            initParams.boundaryObject = inletBoundary;
            nonZeroVFixedDensityILet = new lb::collisions::NonZeroVelocityEquilibriumFixedDensity<lb::kernels::LBGK<
                lb::lattices::D3Q15> >(initParams);

            initParams.boundaryObject = outletBoundary;
            zeroVFixedDensityOLet = new lb::collisions::ZeroVelocityEquilibriumFixedDensity<lb::kernels::LBGK<
                lb::lattices::D3Q15> >(initParams);
            zeroVEqm = new lb::collisions::ZeroVelocityEquilibrium<lb::kernels::LBGK<lb::lattices::D3Q15> >(initParams);
            normal = new lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> >(initParams);
          }

          void tearDown()
          {


            delete inletBoundary;
            delete outletBoundary;

            delete lbgk;

            delete nonZeroVFixedDensityILet;
            delete zeroVFixedDensityOLet;
            delete zeroVEqm;
            delete normal;
            FourCubeBasedTestFixture::tearDown();
          }

          void TestNonZeroVelocityEquilibriumFixedDensity()
          {
            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            geometry::Site dummySite(0, *latDat);
            nonZeroVFixedDensityILet->CalculatePreCollision(hydroVars, dummySite);

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho = inletBoundary->GetBoundaryDensity(0);
            distribn_t expectedV[3];

            LbTestsHelper::CalculateVelocity<lb::lattices::D3Q15>(fOld, expectedV);
            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedV[0],
                                                                  expectedV[1],
                                                                  expectedV[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedV[0],
                                         expectedV[1],
                                         expectedV[2],
                                         expectedFeq,
                                         "Non-0 velocity eqm fixed density, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            nonZeroVFixedDensityILet->Collide(lbmParams, hydroVars);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Non-0 velocity eqm fixed density, collide",
                                                   expectedFeq[ii],
                                                   hydroVars.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

          void TestZeroVelocityEquilibriumFixedDensity()
          {
            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            zeroVFixedDensityOLet->CalculatePreCollision(hydroVars, latDat->GetSite(0));

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho = outletBoundary->GetBoundaryDensity(0);
            distribn_t expectedV[3] = { 0., 0., 0. };

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedV[0],
                                                                  expectedV[1],
                                                                  expectedV[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedV[0],
                                         expectedV[1],
                                         expectedV[2],
                                         expectedFeq,
                                         "0 velocity eqm fixed density, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            zeroVFixedDensityOLet->Collide(lbmParams, hydroVars);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("0 velocity eqm fixed density, collide",
                                                   expectedFeq[ii],
                                                   hydroVars.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

          void TestZeroVelocityEquilibrium()
          {
            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            zeroVEqm->CalculatePreCollision(hydroVars, latDat->GetSite(0));

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho = 0.0;
            distribn_t expectedV[3] = { 0., 0., 0. };

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              expectedRho += fOld[ii];
            }

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedV[0],
                                                                  expectedV[1],
                                                                  expectedV[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedV[0],
                                         expectedV[1],
                                         expectedV[2],
                                         expectedFeq,
                                         "0 velocity eqm, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            zeroVEqm->Collide(lbmParams, hydroVars);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("0 velocity eqm, collide",
                                                   expectedFeq[ii],
                                                   hydroVars.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

          void TestNormal()
          {
            distribn_t allowedError = 1e-10;

            // Initialise the fOld and the hydro vars.
            distribn_t fOld[lb::lattices::D3Q15::NUMVECTORS];

            LbTestsHelper::InitialiseAnisotropicTestData<lb::lattices::D3Q15>(0, fOld);

            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVars(fOld);

            // Test the pre-collision step, which should calculate the correct
            // post-collisional density, velocity and equilibrium distribution.
            normal->CalculatePreCollision(hydroVars, latDat->GetSite(0));

            // Calculate the expected density, velocity and f_eq.
            distribn_t expectedRho;
            distribn_t expectedV[3];

            LbTestsHelper::CalculateRhoVelocity<lb::lattices::D3Q15>(fOld, expectedRho, expectedV);

            distribn_t expectedFeq[lb::lattices::D3Q15::NUMVECTORS];
            LbTestsHelper::CalculateLBGKEqmF<lb::lattices::D3Q15>(expectedRho,
                                                                  expectedV[0],
                                                                  expectedV[1],
                                                                  expectedV[2],
                                                                  expectedFeq);

            // Compare.
            LbTestsHelper::CompareHydros(expectedRho,
                                         expectedV[0],
                                         expectedV[1],
                                         expectedV[2],
                                         expectedFeq,
                                         "Normal, calculate pre collision",
                                         hydroVars,
                                         allowedError);

            // Next, compare the collision function itself. The result should be the equilibrium
            // distribution.
            // Make a copy for the second collision to compare against.
            lb::kernels::HydroVars<lb::kernels::LBGK<lb::lattices::D3Q15> > hydroVarsCopy(hydroVars);

            lbgk->Collide(lbmParams, hydroVars);
            normal->Collide(lbmParams, hydroVarsCopy);

            for (unsigned int ii = 0; ii < lb::lattices::D3Q15::NUMVECTORS; ++ii)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Normal, collide",
                                                   hydroVars.GetFPostCollision()[ii],
                                                   hydroVarsCopy.GetFPostCollision()[ii],
                                                   allowedError);
            }
          }

        private:

          lb::boundaries::BoundaryValues* inletBoundary;
          lb::boundaries::BoundaryValues* outletBoundary;

          lb::kernels::LBGK<lb::lattices::D3Q15>* lbgk;

          lb::collisions::NonZeroVelocityEquilibriumFixedDensity<lb::kernels::LBGK<lb::lattices::D3Q15> >
              * nonZeroVFixedDensityILet;
          lb::collisions::ZeroVelocityEquilibriumFixedDensity<lb::kernels::LBGK<lb::lattices::D3Q15> >
              * zeroVFixedDensityOLet;
          lb::collisions::ZeroVelocityEquilibrium<lb::kernels::LBGK<lb::lattices::D3Q15> >* zeroVEqm;
          lb::collisions::Normal<lb::kernels::LBGK<lb::lattices::D3Q15> >* normal;
      };
      CPPUNIT_TEST_SUITE_REGISTRATION( CollisionTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_COLLISIONTESTS_H_ */
