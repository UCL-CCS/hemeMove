#ifndef HEMELB_VIS_CLUSTERRAYTRACER_H
#define HEMELB_VIS_CLUSTERRAYTRACER_H

//#include <map>
//#include <stack>
//#include <vector>

//#include "constants.h"
#include "geometry/LatticeData.h"
//#include "topology/NetworkTopology.h"

//#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
//#include "vis/VisSettings.h"
#include "vis/Vector3D.h"
#include "vis/XYCoordinates.h"
#include "vis/rayTracer/Cluster.h"
#include "vis/rayTracer/Ray.h"
//#include "vis/rayTracer/ClusterBuilder.h"
//#include "vis/rayTracer/SiteData.h"

namespace hemelb
{
  namespace vis
  {
    namespace raytracer 
    {
      namespace Direction 
      {
	enum Direction
	{
	  X,
	  Y,
	  Z
	};
      }

      class ClusterRayTracer
      {
      public:
	ClusterRayTracer(const Viewpoint& iViewpoint,
			 Screen& iScreen,
			 const DomainStats& iDomainStats,
			 const VisSettings& iVisSettings,
			 const hemelb::geometry::LatticeData& iLatticeData);
	  
	void RenderCluster(const Cluster& iCluster);
      
      private:
	void CalculateSubImage(const Cluster& iCluster);
	
	void UpdateSubImageExtentForCorner
	  (const Vector3D<float>& iCorner,
	   XYCoordinates<float>& ioSubImageLowerLeft,
	   XYCoordinates<float>& ioSubImageUpperRight);

	bool SubImageOffScreen();

	void CropSubImageToScreen();
	
	void CalculateVectorsToClusterSpanAndLowerLeftPixel(const Cluster& iCluster);

	void CastRaysForEachPixel(const Cluster& iCluster);

	void CastRayForPixel(const Cluster& iCluster,
			     const XYCoordinates<int>& iPixel,
			     Vector3D<float> iRayDirection);

	void GetRayUnitsFromViewpointToCluster
	  (const Ray & iRay, float & oMaximumRayUnits,
	   float & oMinimumRayUnits);

	void TraverseVoxels(const Vector3D<float>& block_min,
			    const Vector3D<float>& block_x,
			    const SiteData_t* iSiteData,
			    float t,
			    Ray& bCurrentRay);

	Vector3D<float> CalculateRayUnitsBeforeNextVoxel
	  (const Vector3D<float>& iFirstRayClusterIntersectionToBlockLowerSite,
	   const Vector3D<site_t>& iTruncatedLocationInBlock, const Ray& iRay);

	Vector3D<site_t> RoundToNearestVoxel(const Vector3D<float>& iUnboundLocation);
      
	Direction::Direction DirectionOfLeastTravel
	  (Vector3D<float> iRayUnitsBeforeNextVoxelOrBlock);

	void TraverseBlocks(const Cluster& iCluster,
			    const Vector3D<float>& iLowerSiteToFirstRayClusterIntersection,
			    Ray& iRay);
	
	Vector3D<unsigned int> GetBlockCoordinatesOfFirstIntersectionBlock(
	const Cluster& iCluster,
	Vector3D<float> iLowerSiteToFirstRayClusterIntersection);

	Vector3D<float> CalculateRayUnitsBeforeNextBlock
	  (const Vector3D<float>& lFirstIntersectionToBlockLowerSite,
	   const Ray& iRay);
	
	void UpdateRayData(const SiteData_t* iSiteData,
			   float ray_t,
			   float ray_segment,
			   Ray* bCurrentRay);

	void UpdateColour(float dt, const float palette[3], float col[3]);

	const Viewpoint& mViewpoint;

	Screen& mScreen;

	const DomainStats& mDomainStats;

	const VisSettings& mVisSettings;

	const hemelb::geometry::LatticeData& mLatticeData;

	Vector3D<float> mCameraToBottomLeftPixel;

	Vector3D<float> mLowerSiteCordinatesOfClusterRelativeToViewpoint;

	XYCoordinates<int> mSubImageLowerLeftPixelCoordinates;
	
	XYCoordinates<int> mSubImageUpperRightPixelCoordinates;

	//Vectors from the viewpoint centre 
	//to the maximum and minimum site span
	//locations respectively
	//(Formerly AABB)
	Vector3D<float> mViewpointCentreToMaxSite;
	Vector3D<float> mViewpointCentreToMinSite;

	
      };
    }
  }
}

#endif // HEMELB_VIS_CLUSTERRENDERER_H
