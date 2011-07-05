#ifndef HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGK_H
#define HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGK_H

#include "lb/collisions/implementations/CollisionOperator.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      namespace implementations
      {

        class LBGK : public CollisionOperator
        {

          public:
            static void getSiteValues(const distribn_t* f,
                                      distribn_t &density,
                                      distribn_t &v_x,
                                      distribn_t &v_y,
                                      distribn_t &v_z,
                                      distribn_t* f_eq);

            static void getBoundarySiteValues(const distribn_t* f,
                                              const distribn_t &density,
                                              const distribn_t &v_x,
                                              const distribn_t &v_y,
                                              const distribn_t &v_z,
                                              distribn_t* f_eq);

            static distribn_t getOperatorElement(distribn_t &f_i,
                                                 distribn_t &f_eq_i,
                                                 const LbmParameters* iLbmParams);

        };

      }
    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGK_H */
