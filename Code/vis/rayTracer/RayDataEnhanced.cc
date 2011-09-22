#include "vis/rayTracer/RayDataEnhanced.h"

namespace hemelb
{
  namespace vis
  { 
    namespace raytracer 
    {
      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mSurfaceNormalLightnessRange = 0.3F;
     
      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mParallelSurfaceAttenuation = 0.5F;
      
      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mLowestLightness = 0.3F;

      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mVelocityHueMin = 240.0F;
     
      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mVelocityHueRange = 120.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mVelocitySaturation = 1.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mStressHue = 230.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mStressSaturationRange = 0.5F;

      template <>
      const float RayDataEnhanced<DepthCuing::MIST>::mStressSaturationMin = 0.5F;    
     
      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mSurfaceNormalLightnessRange = 0.5F;

      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mParallelSurfaceAttenuation = 0.75F;
      
      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mLowestLightness = 0.3F;
      
      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mVelocityHueMin = 240.0F;
     
      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mVelocityHueRange = 120.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mVelocitySaturation = 1.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mStressHue = 230.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mStressSaturationRange = 0.5F;

      template <>
      const float RayDataEnhanced<DepthCuing::NONE>::mStressSaturationMin = 0.5F;    

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mSurfaceNormalLightnessRange = 0.3F;

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mParallelSurfaceAttenuation = 0.5F;
      
      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mLowestLightness = 0.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mVelocityHueMin = 240.0F;
     
      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mVelocityHueRange = 120.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mVelocitySaturation = 1.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mStressHue = 230.0F;

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mStressSaturationRange = 0.5F;

      template <>
      const float RayDataEnhanced<DepthCuing::DARKNESS>::mStressSaturationMin = 0.5F;    
    }
  }

  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced
  <vis::raytracer::DepthCuing::MIST> >
  ::RegisterMpiDataType()
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
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced
				 <vis::raytracer::DepthCuing::DARKNESS> >::RegisterMpiDataType()
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
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::raytracer::RayDataEnhanced
				 <vis::raytracer::DepthCuing::NONE> >::RegisterMpiDataType()
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
