#ifndef HEMELB_VIS_RAYTRACER_H
#define HEMELB_VIS_RAYTRACER_H

#include "constants.h"
#include "geometry/LatticeData.h"
#include "topology/NetworkTopology.h"

#include "vis/DomainStats.h"
#include "vis/Screen.h"
#include "vis/Viewpoint.h"
#include "vis/VisSettings.h"

namespace hemelb
{
  namespace vis
  {
    struct Location
    {
      site_t i, j, k;
      
      Location() {};

      Location(site_t iI, site_t iJ, site_t iK)
      {
	i = iI;
	j = iJ;
	k = iK;
      }

      Location(site_t iX)
      {
	i=iX;
	j=iX;
	k=iX;
      }

      Location operator+(const Location right)
      {
	return Location(this->i + right.i,
			this->j + right.j,
			this->k + right.k);
      }

      Location operator*(const site_t multiplier)
      {
	return Location(this->i * multiplier,
			this->j * multiplier,
			this->k * multiplier);
      }

    };

    struct Cluster
    {
      float minmax_x[2], minmax_y[2], minmax_z[2];
      
      //Stores the lowest x, y and z point of the Cluster 
      float blockCoordinates[3];
      
      //Stores the size of the cluster
      unsigned short int blocks_x;
      unsigned short int blocks_y;
      unsigned short int blocks_z;
    };

    class RayTracer
    {
      public:
        // Constructor and destructor do all the usual stuff.
        RayTracer(const geometry::LatticeData* iLatDat,
                  const DomainStats* iDomainStats,
                  Screen* iScreen,
                  Viewpoint* iViewpoint,
                  VisSettings* iVisSettings);
        ~RayTracer();

        // Method to update the voxel corresponding to site i with its
        // newly calculated density, velocity and stress.
        void UpdateClusterVoxel(site_t i,
                                distribn_t density,
                                distribn_t velocity,
                                distribn_t stress);

        // Render the current state into an image.
        void Render();

      private:

        struct Ray
        {
            float Direction[3];
            float InverseDirection[3];
            float Length;

            float VelocityColour[3];
            float StressColour[3];
            float Stress;
            float Density;
            float MinT;
        };

        struct AABB
        {
            float acc_1, acc_2, acc_3, acc_4, acc_5, acc_6;
        };


        void UpdateRayData(const float flow_field[3],
                           float ray_t,
                           float ray_segment,
                           Ray* bCurrentRay);

        void TraverseVoxels(const float block_min[3],
                            const float block_x[3],
                            const float voxel_flow_field[],
                            float t,
                            Ray* bCurrentRay,
                            const bool xyz_is_1[3]);

        void TraverseBlocks(const Cluster* cluster,
                            const bool xyz_Is_1[3],
                            const float ray_dx[3],
                            float **block_flow_field,
                            Ray *bCurrentRay);

        void AABBvsRay(const AABB* aabb,
                       const float inverseDirection[3],
                       const bool xyzComponentIsPositive[3],
                       float* t_near,
                       float* t_far);

        void UpdateColour(float dt, const float palette[3], float col[3]);

        void BuildClusters();

        const geometry::LatticeData* mLatDat;

        const DomainStats* mDomainStats;
        Screen* mScreen;
        Viewpoint* mViewpoint;
        VisSettings* mVisSettings;

        std::vector<Cluster> mClusters;
        float **cluster_voxel;
        float ***cluster_flow_field;

        const float mBlockSizeFloat;
        const float mBlockSizeInverse;
        const site_t block_size2, block_size3, block_size_1;
        const site_t blocks_yz;
    };

  }
}

#endif // HEMELB_VIS_RAYTRACER_H
