#ifndef HEMELB_LB_LOCALLATTICEDATA_H
#define HEMELB_LB_LOCALLATTICEDATA_H

#include "constants.h"
#include "D3Q15.h"
#include "debug/Debugger.h"

namespace hemelb
{
  namespace lb
  {
    enum SiteType
    {
      SOLID_TYPE = 0U,
      FLUID_TYPE = 1U,
      INLET_TYPE = 2U,
      OUTLET_TYPE = 3U
    };

    class LocalLatticeData
    {
      public:
        LocalLatticeData(int iLocalFluidSites);
        ~LocalLatticeData();

        int GetStreamedIndex(int iSiteIndex, int iDirectionIndex) const;
        double GetCutDistance(int iSiteIndex, int iDirection) const;
        bool HasBoundary(int iSiteIndex, int iDirection) const;
        int GetBoundaryId(int iSiteIndex) const;
        const double *GetNormalToWall(int iSiteIndex) const;
        SiteType GetSiteType(int iSiteIndex) const;
        int GetLocalFluidSiteCount() const;

        void SetNeighbourLocation(unsigned int iSiteIndex,
                                  unsigned int iDirection,
                                  unsigned int iValue);
        void SetWallNormal(int iSiteIndex, const double iNormal[3]);
        void SetDistanceToWall(int iSiteIndex, const double iCutDistance[D3Q15::NUMVECTORS - 1]);

        void SetSharedSiteCount(int iSharedCount);
      public:
        int my_inner_sites;
        int my_inner_collisions[COLLISION_TYPES];
        int my_inter_collisions[COLLISION_TYPES];

        double *FOld;
        double *FNew;

        // TODO sadly this has to be public, due to some budgetry in the way we determine site type.
        // SiteType || FluidSite and SiteType && FluidSite have different significances...
        unsigned int *mSiteData;

      private:
        int LocalFluidSites;
        unsigned int *mFNeighbours;
        double *mDistanceToWall;
        double *mWallNormalAtSite;
    };
  }
}

#endif /* HEMELB_LB_LOCALLATTICEDATA_H */
