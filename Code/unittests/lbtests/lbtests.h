#ifndef HEMELB_UNITTESTS_LBTESTS_LBTESTS_H
#define HEMELB_UNITTESTS_LBTESTS_LBTESTS_H

#include <cppunit/TestCaller.h>
#include <cppunit/TestFixture.h>

#include "unittests/lbtests/KernelTests.h"
#include "unittests/lbtests/CollisionTests.h"
#include "unittests/lbtests/StreamerTests.h"

namespace hemelb
{
  namespace unittests
  {
    namespace lbtests
    {

      /**
       * Class containing tests for the functionality of the lattice-Boltzmann kernels.
       */
      class LbTestSuite : public CppUnit::TestSuite
      {
        public:
          LbTestSuite()
          {
            addTest(new CppUnit::TestCaller<KernelTests>("TestEntropicCalculationsAndCollision",
                                                         &KernelTests::TestEntropicCalculationsAndCollision));
            addTest(new CppUnit::TestCaller<KernelTests>("TestLBGKCalculationsAndCollision",
                                                         &KernelTests::TestLBGKCalculationsAndCollision));

            addTest(new CppUnit::TestCaller<CollisionTests>("TestNonZeroVelocityEquilibriumFixedDensity",
                                                            &CollisionTests::TestNonZeroVelocityEquilibriumFixedDensity));
            addTest(new CppUnit::TestCaller<CollisionTests>("TestZeroVelocityEquilibriumFixedDensity",
                                                            &CollisionTests::TestZeroVelocityEquilibriumFixedDensity));
            addTest(new CppUnit::TestCaller<CollisionTests>("TestZeroVelocityEquilibrium",
                                                            &CollisionTests::TestZeroVelocityEquilibrium));
            addTest(new CppUnit::TestCaller<CollisionTests>("TestNormal",
                                                            &CollisionTests::TestNormal));

            addTest(new CppUnit::TestCaller<StreamerTests>("TestSimpleCollideAndStream",
                                                           &StreamerTests::TestSimpleCollideAndStream));
          }
      };

    }
  }
}

#endif /* HEMELB_UNITTESTS_LBTESTS_LBTESTS_H */
