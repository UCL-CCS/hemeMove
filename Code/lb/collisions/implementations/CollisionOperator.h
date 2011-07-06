#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_COLLISIONOPERATOR_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_COLLISIONOPERATOR_H

#include "lb/collisions/HFunction.h"
#include "util/utilityFunctions.h"
#include "constants.h"
#include "lb/LbmParameters.h"
#include "D3Q15.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        class CollisionOperator
        {

          protected:
            CollisionOperator();

            static double getAlpha(const distribn_t* lFOld,const  distribn_t* lFEq, double prevAlpha);

        };

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_COLLISIONOPERATOR_H */
