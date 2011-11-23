#ifndef HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H
#define HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H

#include "util/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      // Class to hold information about the velocity field at some point.
      class VelocitySiteData
      {
        public:
          VelocitySiteData()
          {
            proc_id = -1;
            counter = -1;
          }

          /**
           * Counter describes how many iterations have passed since the objects creation.
           * It allows for quickly checking whether the local velocity field has been
           * calculated this iteration.
           **/
          site_t counter;

          /**
           * The rank this volume of the geometry lives on.
           */
          proc_t proc_id;

          /**
           * The site id of the volume represented.
           */
          site_t site_id;

          /**
           * The velocity calculated for that site.
           */
          util::Vector3D<float> velocity;
      };
    }
  }
}

#endif //HEMELB_VIS_STREAKLINEDRAWER_VELOCITYSITEDATA_H
