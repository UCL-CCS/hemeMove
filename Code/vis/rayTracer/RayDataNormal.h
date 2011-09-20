#ifndef HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
#define HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H

#include "mpiInclude.h"
#include "vis/DomainStats.h"
#include "vis/rayTracer/RayData.h"
#include "vis/VisSettings.h"
#include "vis/Vector3D.h"
#include "lb/LbmParameters.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer
    { 
      //NB functions prefixed Do should only be called by the base class
      class RayDataNormal : public RayData<RayDataNormal>
      {
      public:
	RayDataNormal();
	
        //Used to process the ray data for a normal (non-wall) fluid site
	void DoUpdateDataForNormalFluidSite(const SiteData_t& iSiteData, 
					    const Vector3D<float>& iRayDirection,
					    const float iRayLengthInVoxel,
					    const DomainStats& iDomainStats,
					    const VisSettings& iVisSettings);
	//Passed as references since pointer can't be
	//meaningly transfered over MPI

	//Used to process the ray data for wall site
	void DoUpdateDataForWallSite(const SiteData_t& iSiteData, 
				     const Vector3D<float>& iRayDirection,
				     const float iRayLengthInVoxel,
				     const DomainStats& iDomainStats,
				     const VisSettings& iVisSettings,
				     const double* iWallNormal);
	
	//Carries out the merging of the ray data in this
	//inherited type, for different segments of the same ray
	void DoMergeIn(const RayDataNormal& iOtherRayData,
		       const VisSettings& iVisSettings);

	//Obtains the colour representing the velocity ray trace
	void DoGetVelocityColour(unsigned char oColour[3]) const;
 
	//Obtains the colour represting the stress ray trace
	void DoGetStressColour(unsigned char oColour[3]) const;

      private:
	void UpdateVelocityColour(float iDt, const float iPalette[3]);

	void UpdateStressColour(float iDt, const float iPalette[3]);
      
	float mVelR, mVelG, mVelB;
	float mStressR, mStressG, mStressB;
      };	
    }
  }

}

#endif // HEMELB_VIS_RAYTRACER_RAYDATANORMAL_H
