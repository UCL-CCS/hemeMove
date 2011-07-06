#include "lb/collisions/implementations/CollisionOperator.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        CollisionOperator::CollisionOperator()
        {

        }

        double CollisionOperator::getAlpha(const distribn_t* lF,
                                           const distribn_t* lFEq,
                                           double prevAlpha)
        {
          bool small = true;

          for (unsigned int i = 0; i < D3Q15::NUMVECTORS; i++)
          {
            // Papers suggest f_eq - f < 0.001 or (f_eq - f)/f < 0.01 for the point to have approx alpha = 2
            // Accuracy can change depending on stability requirements, because the more NR evaluations it skips
            // the more of the simulation is in the LBGK limit.
            if (fabs(lFEq[i] - lF[i]) > 1.0E-2)
            {
              small = false;
              break;
            }
          }

          if (small)
            return 2.0;

          HFunction HFunc(lF, lFEq);

          // Accuracy is set to 1.0E-3 as this works for difftest.
          return (hemelb::util::NumericalFunctions::NewtonRaphson(&HFunc, prevAlpha, 1.0E-3));

        }

      }
    }
  }
}

