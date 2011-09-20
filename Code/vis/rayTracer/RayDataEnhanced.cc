#include "vis/rayTracer/RayDataEnhanced.h"

namespace hemelb
{
  namespace vis
  { 
    namespace raytracer 
    {
      template <>
      const float RayDataEnhanced<true>::mSurfaceNormalLightnessRange = 0.3F;
     
      template <>
      const float RayDataEnhanced<true>::mParallelSurfaceAttenuation = 0.5F;
      
      template <>
      const float RayDataEnhanced<true>::mLowestLightness = 0.3F;

      template <>
      const float RayDataEnhanced<false>::mSurfaceNormalLightnessRange = 0.5F;

      template <>
      const float RayDataEnhanced<false>::mParallelSurfaceAttenuation = 0.75F;
      
      template <>
      const float RayDataEnhanced<false>::mLowestLightness = 0.3F;
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced<false> >::RegisterMpiDataType()
  {
    int lRayDataEnhancedCount = 8;
    int lRayDataEnhancedBlocklengths[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype lRayDataEnhancedTypes[8] = { MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_UB };

    MPI_Aint lRayDataEnhancedDisps[8];

    lRayDataEnhancedDisps[0] = 0;

    for (int i = 1; i < lRayDataEnhancedCount; i++)
    {
      if (lRayDataEnhancedTypes[i - 1] == MPI_FLOAT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1]
	  + (sizeof(float) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_INT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(int) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_UNSIGNED)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(unsigned) * lRayDataEnhancedBlocklengths[i - 1]);
      }
    }
    MPI_Datatype type;
    MPI_Type_struct(lRayDataEnhancedCount,
		    lRayDataEnhancedBlocklengths,
		    lRayDataEnhancedDisps,
		    lRayDataEnhancedTypes,
		    &type);
    MPI_Type_commit(&type);
    return type;
  }

template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced<true> >::RegisterMpiDataType()
  {
    int lRayDataEnhancedCount = 8;
    int lRayDataEnhancedBlocklengths[8] = { 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype lRayDataEnhancedTypes[8] = { MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_FLOAT,
					      MPI_UB };

    MPI_Aint lRayDataEnhancedDisps[8];

    lRayDataEnhancedDisps[0] = 0;

    for (int i = 1; i < lRayDataEnhancedCount; i++)
    {
      if (lRayDataEnhancedTypes[i - 1] == MPI_FLOAT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1]
	  + (sizeof(float) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_INT)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(int) * lRayDataEnhancedBlocklengths[i - 1]);
      }
      else if (lRayDataEnhancedTypes[i - 1] == MPI_UNSIGNED)
      {
	lRayDataEnhancedDisps[i] = lRayDataEnhancedDisps[i - 1] + (sizeof(unsigned) * lRayDataEnhancedBlocklengths[i - 1]);
      }
    }
    MPI_Datatype type;
    MPI_Type_struct(lRayDataEnhancedCount,
		    lRayDataEnhancedBlocklengths,
		    lRayDataEnhancedDisps,
		    lRayDataEnhancedTypes,
		    &type);
    MPI_Type_commit(&type);
    return type;
  }
}
