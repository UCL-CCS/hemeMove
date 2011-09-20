//#define NDEBUG;
#include <assert.h>

#include "geometry/LatticeData.h"
#include "vis/Vector3D.h"
#include "vis/BlockTraverser.h"
#include "vis/rayTracer/RayTracer.h"

namespace hemelb
{
  namespace vis
  {
    BlockTraverser::BlockTraverser
    (const geometry::LatticeData& iLatDat)
      : VolumeTraverser(),
	mLatticeData(iLatDat)
    { 		
    }

    BlockTraverser::~BlockTraverser()
    {
    }
      
    site_t BlockTraverser::CurrentBlockNumber()
    {
      return GetCurrentIndex();
    }

    Vector3D<site_t> BlockTraverser::GetSiteCoordinatesOfLowestSiteInCurrentBlock()
    {
      return GetCurrentLocation()*mLatticeData.GetBlockSize();
    }
      		
   	    
    geometry::LatticeData::BlockData *
    BlockTraverser::GetCurrentBlockData()
    {
      return mLatticeData.GetBlock(mCurrentNumber);
    }

    geometry::LatticeData::BlockData *
    BlockTraverser::GetBlockDataForLocation
    (const Vector3D<site_t>& iLocation)
    {
      return mLatticeData.GetBlock(GetIndexFromLocation(iLocation));
    }      

    site_t BlockTraverser::GetBlockSize()
    {
      return mLatticeData.GetBlockSize();
    }

    SiteTraverser 
    BlockTraverser::GetSiteTraverser()
    {
      return SiteTraverser(mLatticeData);
    }	      

    bool BlockTraverser::IsValidLocation(Vector3D<site_t> iBlock)
    {
      return mLatticeData.IsValidBlockSite
	(iBlock.x,
	 iBlock.y,
	 iBlock.z);
    }
  
    bool BlockTraverser::GoToNextBlock()
    {
      return TraverseOne();
    }
	  
    site_t BlockTraverser::GetXCount() 
    {
      return mLatticeData.GetXBlockCount();
    }

    site_t BlockTraverser::GetYCount()
    {
      return mLatticeData.GetYBlockCount();
    }

    site_t BlockTraverser::GetZCount() 
    {
      return mLatticeData.GetZBlockCount();
    }
  }
}
