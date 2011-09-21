#ifndef HEMELB_VIS_RAYTRACER_VOLUMETRAVERSER_H
#define HEMELB_VIS_RAYTRACER_VOLUMETRAVERSER_H

#include "util/Vector3D.h"

namespace hemelb
{
  namespace vis
  {
    //Volume traverse is used to sequentially traverse a 
    //3D structure maintaining the index and Location 
    //within volume
    class VolumeTraverser
    {
    public:
      VolumeTraverser();

      util::Vector3D<site_t> GetCurrentLocation();
			
      void SetCurrentLocation(const util::Vector3D<site_t>& iLocation);

      site_t GetX();
      
      site_t GetY();
      
      site_t GetZ();
      
      site_t GetCurrentIndex();

      site_t GetIndexFromLocation(util::Vector3D<site_t> iLocation);

      //Increments the index by one and update the location accordingly
      //Returns true if successful or false if the whole volume has been
      //traversed
      bool TraverseOne();

      void IncrementX();
      void IncrementY();
      void IncrementZ();
	
      void DecrementX();
      void DecrementY();
      void DecrementZ();

      bool CurrentLocationValid();

      //Virtual methods which must be defined for correct traversal
      virtual site_t GetXCount() = 0;
      virtual site_t GetYCount() = 0;
      virtual site_t GetZCount() = 0;
	    
    protected:
      util::Vector3D<site_t> mCurrentLocation;
      site_t mCurrentNumber;
    };
  }
}


#endif // HEMELB_VIS_RAYTRACER_VOLUMETRAVERSER_H
