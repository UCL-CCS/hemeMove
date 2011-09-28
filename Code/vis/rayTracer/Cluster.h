#ifndef HEMELB_VIS_RAYTRACER_CLUSTER_H
#define HEMELB_VIS_RAYTRACER_CLUSTER_H

#include <vector>

#include "util/Vector3D.h"

#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    {
      //The cluster structure stores data relating to the clusters
      //used by the RayTracer, in an optimal format
      //Cluster are produced by the ClusterSharedFactory
      //Caution: the data within the flow field is altered by means
      //of pointers obtained from the GetClusterSharedVoxelDataPointer
      //method
      template<typename Derived>
      class Cluster
      {
        public:
          unsigned int GetBlockIdFrom3DBlockLocation(util::Vector3D<unsigned int> iLocation) const
          {
            return ((Derived*)(this))->DoGetBlockIdFrom3DBlockLocation(iLocation);
          }

          //Resizes the vectors so as to be the correct size based on the stored sizes
          void ResizeVectors()
          {
            ((Derived*)(this))->DoResizeVectors();
          }

          void ResizeVectorsForBlock(site_t iBlockNumber, site_t iSize)
          {
            ((Derived*)(this))->DoResizeVectorsForBlock(iBlockNumber, iSize);
          }

          //Returns true if there is site data for a given block
          bool BlockContainsSites(site_t iBlockNumber) const
          {
            return ((const Derived*)(this))->DoBlockContainsSites(iBlockNumber);
          }

          //Get SiteData arary for site
          const SiteData_t* GetSiteData(site_t iBlockNumber) const
          {
            return ((const Derived*)(this))->DoGetSiteData(iBlockNumber);
          }

          const SiteData_t* GetSiteData(site_t iBlockNumber, site_t iSiteNumber) const
          {
            return ((const Derived*)(this))->DoGetSiteData(iBlockNumber, iSiteNumber);
          }

          const double* GetWallData(site_t iBlockNumber, site_t iSiteNumber) const
          {
            return ((const Derived*)(this))->DoGetWallData(iBlockNumber, iSiteNumber);
          }

          void SetWallData(site_t iBlockNumber, site_t iSiteNumber, double* iData)
          {
            return ((Derived*)(this))->DoSetWallData(iBlockNumber, iSiteNumber, iData);
          }

          static bool NeedsWallNormals()
          {
            return Derived::DoNeedsWallNormals();
          }
      };
    }
  }
}

#endif // HEMELB_VIS_RAYTRACER_CLUSTER_H
