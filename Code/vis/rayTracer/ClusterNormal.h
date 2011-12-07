#ifndef HEMELB_VIS_RAYTRACER_CLUSTERNORMAL_H
#define HEMELB_VIS_RAYTRACER_CLUSTERNORMAL_H

#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      class ClusterNormal : public Cluster<ClusterNormal>
      {
        public:
          ClusterNormal(unsigned short xBlockCount,
                        unsigned short yBlockCount,
                        unsigned short zBlockCount,
                        const util::Vector3D<float>& minimalSite,
                        const util::Vector3D<float>& maximalSite,
                        const util::Vector3D<float>& minimalSiteOnMinimalBlock);

          const double* DoGetWallData(site_t iBlockNumber, site_t iSiteNumber) const;

          void DoSetWallData(site_t iBlockNumber, site_t iSiteNumber, const double* const iData);
      };

    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTERNORMAL_H
