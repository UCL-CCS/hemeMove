#ifndef HEMELB_UNITTESTS_LBTESTS_LATTICETESTS_H
#define HEMELB_UNITTESTS_LBTESTS_LATTICETESTS_H

#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include "lb/lattices/D3Q15.h"
#include "lb/lattices/D3Q19.h"
#include "lb/lattices/D3Q27.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {
      class LatticeTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (LatticeTests);
          CPPUNIT_TEST (TestD3Q15);
          CPPUNIT_TEST (TestD3Q19);
          CPPUNIT_TEST (TestD3Q27);CPPUNIT_TEST_SUITE_END();

        public:

          LatticeTests() :
              epsilon(1e-10)
          {
          }
          // had to move this here for compilation portability, see http://stackoverflow.com/questions/370283/why-cant-i-have-a-non-integral-static-const-member-in-a-class
          // wasn't compiling under cray compilers)

          void setUp()
          {
          }

          void tearDown()
          {

          }

          void TestD3Q15()
          {
            TestLattice<lb::lattices::D3Q15>();
          }

          void TestD3Q19()
          {
            TestLattice<lb::lattices::D3Q19>();
          }

          void TestD3Q27()
          {
            TestLattice<lb::lattices::D3Q27>();
          }

        private:
          template<class LatticeType>
          void TestLattice()
          {
            /*
             Interface being tested:

             static const unsigned int NUMVECTORS
             static const int CX[NUMVECTORS];
             static const int CY[NUMVECTORS];
             static const int CZ[NUMVECTORS];

             Require that each vector element is in {0,1,-1} and that each vector is unique.
             */

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              CPPUNIT_ASSERT(LatticeType::CX[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CY[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CZ[direction] >= -1);
              CPPUNIT_ASSERT(LatticeType::CX[direction] <= 1);
              CPPUNIT_ASSERT(LatticeType::CY[direction] <= 1);
              CPPUNIT_ASSERT(LatticeType::CZ[direction] <= 1);

              for (Direction otherDirection = 0; otherDirection < LatticeType::NUMVECTORS; ++otherDirection)
              {
                if (otherDirection == direction)
                {
                  continue;
                }

                CPPUNIT_ASSERT(LatticeType::CX[direction] != LatticeType::CX[otherDirection] || LatticeType::CY[direction] != LatticeType::CY[otherDirection] || LatticeType::CZ[direction] != LatticeType::CZ[otherDirection]);
              }
            }

            /*
             static const int INVERSEDIRECTIONS[NUMVECTORS];

             Require that inverses are in pairs (this inherently checks for uniqueness, obvs)
             */
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
              Direction inverse = LatticeType::INVERSEDIRECTIONS[direction];

              CPPUNIT_ASSERT(direction == LatticeType::INVERSEDIRECTIONS[inverse]);
            }

            /*
             static void CalculateDensityAndVelocity(const distribn_t f[],
             distribn_t &density,
             distribn_t &v_x,
             distribn_t &v_y,
             distribn_t &v_z);
             */

            distribn_t f_data[LatticeType::NUMVECTORS];

            // The 3 here is essentially a seed that relates to the magnitude of the density.
            LbTestsHelper::InitialiseAnisotropicTestData<LatticeType>(3, f_data);

            distribn_t density, velocity[3], expectedDensity, expectedVelocity[3];
            LatticeType::CalculateDensityAndVelocity(f_data, density, velocity[0], velocity[1], velocity[2]);

            LbTestsHelper::CalculateRhoVelocity<LatticeType>(f_data, expectedDensity, expectedVelocity);

            CPPUNIT_ASSERT(density == expectedDensity);

            CPPUNIT_ASSERT(velocity[0] == expectedVelocity[0]);
            CPPUNIT_ASSERT(velocity[1] == expectedVelocity[1]);
            CPPUNIT_ASSERT(velocity[2] == expectedVelocity[2]);

            /*
             static void CalculateFeq(const distribn_t &density,
             const distribn_t &xMomentum,
             const distribn_t &yMomentum,
             const distribn_t &zMomentum,
             distribn_t f_eq[]);

             static void CalculateEntropicFeq(const distribn_t &density,
             const distribn_t &v_x,
             const distribn_t &v_y,
             const distribn_t &v_z,
             distribn_t f_eq[]);
             */
            distribn_t equilibriumF[LatticeType::NUMVECTORS], expectedEquilibriumF[LatticeType::NUMVECTORS],
                equilibriumEntropicFAnsumali[LatticeType::NUMVECTORS],
                expectedEquilibriumEntropicFAnsumali[LatticeType::NUMVECTORS],
                equilibriumEntropicFChikatamarla[LatticeType::NUMVECTORS];

            // These values chosen as they're pairwise coprime. Probably doesn't matter.
            distribn_t targetDensity = 0.95, targetH[3] = { 0.002, 0.003, 0.004 };

            LatticeType::CalculateFeq(targetDensity, targetH[0], targetH[1], targetH[2], equilibriumF);
            LatticeType::CalculateEntropicFeqAnsumali(targetDensity,
                                                      targetH[0],
                                                      targetH[1],
                                                      targetH[2],
                                                      equilibriumEntropicFAnsumali);
            LatticeType::CalculateEntropicFeqChik(targetDensity,
                                                  targetH[0],
                                                  targetH[1],
                                                  targetH[2],
                                                  equilibriumEntropicFChikatamarla);

            LbTestsHelper::CalculateLBGKEqmF<LatticeType>(targetDensity,
                                                          targetH[0],
                                                          targetH[1],
                                                          targetH[2],
                                                          expectedEquilibriumF);

            LbTestsHelper::CalculateAnsumaliEntropicEqmF<LatticeType>(targetDensity,
                                                                      targetH[0],
                                                                      targetH[1],
                                                                      targetH[2],
                                                                      expectedEquilibriumEntropicFAnsumali);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedEquilibriumF[direction], equilibriumF[direction], epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedEquilibriumEntropicFAnsumali[direction],
                                           equilibriumEntropicFAnsumali[direction],
                                           epsilon);
            }

            /*
             * It's also the case that these should be invertible (i.e. the density and velocity should be what we started with).
             */
            distribn_t entropicCalculatedDensityAnsumali, entropicCalculatedVelocityAnsumali[3],
                entropicCalculatedDensityChikatamarla, entropicCalculatedVelocityChikatamarla[3], calculatedDensity,
                calculatedVelocity[3];

            LatticeType::CalculateDensityAndVelocity(equilibriumF,
                                                     calculatedDensity,
                                                     calculatedVelocity[0],
                                                     calculatedVelocity[1],
                                                     calculatedVelocity[2]);

            LatticeType::CalculateDensityAndVelocity(equilibriumEntropicFAnsumali,
                                                     entropicCalculatedDensityAnsumali,
                                                     entropicCalculatedVelocityAnsumali[0],
                                                     entropicCalculatedVelocityAnsumali[1],
                                                     entropicCalculatedVelocityAnsumali[2]);

            LatticeType::CalculateDensityAndVelocity(equilibriumEntropicFChikatamarla,
                                                     entropicCalculatedDensityChikatamarla,
                                                     entropicCalculatedVelocityChikatamarla[0],
                                                     entropicCalculatedVelocityChikatamarla[1],
                                                     entropicCalculatedVelocityChikatamarla[2]);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(calculatedDensity, targetDensity, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedDensityAnsumali, targetDensity, epsilon);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedDensityChikatamarla, targetDensity, epsilon);

            for (Direction direction = 0; direction < 3; direction++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(calculatedVelocity[direction], targetH[direction], epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedVelocityAnsumali[direction], targetH[direction], epsilon);
              CPPUNIT_ASSERT_DOUBLES_EQUAL(entropicCalculatedVelocityChikatamarla[direction],
                                           targetH[direction],
                                           epsilon);
            }

            /*
             * static LatticeInfo* GetLatticeInfo();
             *
             * This must have the same values as the corresponding properties on the main lattice object.
             */
            lb::lattices::LatticeInfo& latticeInfo = LatticeType::GetLatticeInfo();

            CPPUNIT_ASSERT(latticeInfo.GetNumVectors() == LatticeType::NUMVECTORS);

            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; direction++)
            {
              CPPUNIT_ASSERT(latticeInfo.GetInverseIndex(direction) == LatticeType::INVERSEDIRECTIONS[direction]);

              const util::Vector3D<int>& velocityVector = latticeInfo.GetVector(direction);

              for (Direction index = 0; index < 3; index++)
              {
                CPPUNIT_ASSERT(velocityVector[index] == LatticeType::discreteVelocityVectors[index][direction]);
              }
            }

            /*
             * TODO: Currently untested
             * * CalculateDensityVelocityFEq (it should just call other functions that *are* tested)
             * * CalculateEntropicDensityVelocityFEq (as above)
             * * CalculateVonMisesStress (probably needs a manually calculated test)
             * * CalculateShearStress (as above)
             * * CalculatePiTensor(as above)
             * * CalculateShearRate (as above)
             *
             */
          }

          const double epsilon;
      };

      CPPUNIT_TEST_SUITE_REGISTRATION (LatticeTests);
    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_LATTICETESTS_H */
