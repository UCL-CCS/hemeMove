#include <math.h>
#include <vector>
#include <cassert>

#include "geometry/BlockTraverser.h"
#include "geometry/SiteTraverser.h"
#include "vis/streaklineDrawer/VelocityField.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      VelocityField::VelocityField(std::map<proc_t, NeighbouringProcessor>& neighbouringProcessorsIn) :
        counter(0), neighbouringProcessors(neighbouringProcessorsIn)
      {
      }

      void VelocityField::BuildVelocityField(const geometry::LatticeData& latDat)
      {
        velocityField.resize(latDat.GetBlockCount());

        // Iterate over each block with some sites on this rank.
        geometry::BlockTraverser blockTraverser(latDat);
        do
        {
          geometry::BlockData* block = blockTraverser.GetCurrentBlockData();

          if (block->site_data == NULL)
          {
            continue;
          }

          geometry::SiteTraverser siteTraverser(latDat);
          do
          {
            // Only interested if the site lives on this rank.
            if (topology::NetworkTopology::Instance()->GetLocalRank()
                != block->ProcessorRankForEachBlockSite[siteTraverser.GetCurrentIndex()])
            {
              continue;
            }

            // Calculate the bounds of the unit cube around the current site (within the lattice)
            const site_t startI =
                util::NumericalFunctions::max<site_t>(0,
                                                      blockTraverser.GetX()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetX() - 1);

            const site_t startJ =
                util::NumericalFunctions::max<site_t>(0,
                                                      blockTraverser.GetY()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetY() - 1);

            const site_t startK =
                util::NumericalFunctions::max<site_t>(0,
                                                      blockTraverser.GetZ()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetZ() - 1);

            const site_t endI =
                util::NumericalFunctions::min<site_t>(latDat.GetXSiteCount() - 1,
                                                      blockTraverser.GetX()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetX() + 1);

            const site_t endJ =
                util::NumericalFunctions::min<site_t>(latDat.GetYSiteCount() - 1,
                                                      blockTraverser.GetY()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetY() + 1);

            const site_t endK =
                util::NumericalFunctions::min<site_t>(latDat.GetZSiteCount() - 1,
                                                      blockTraverser.GetZ()
                                                          * blockTraverser.GetBlockSize()
                                                          + siteTraverser.GetZ() + 1);

            // Iterate over the sites in the unit cube.
            for (site_t neighbourI = startI; neighbourI <= endI; neighbourI++)
            {
              for (site_t neighbourJ = startJ; neighbourJ <= endJ; neighbourJ++)
              {
                for (site_t neighbourK = startK; neighbourK <= endK; neighbourK++)
                {
                  // Get the rank that the neighbour lives on.
                  const proc_t* neigh_proc_id = latDat.GetProcIdFromGlobalCoords(util::Vector3D<
                      site_t>(neighbourI, neighbourJ, neighbourK));

                  // If we have data for it, we should initialise a block in the velocity field
                  // for the neighbour site.
                  if (neigh_proc_id == NULL || *neigh_proc_id == BIG_NUMBER2)
                  {
                    continue;
                  }

                  InitializeVelocityFieldBlock(latDat,
                                               util::Vector3D<site_t>(neighbourI,
                                                                      neighbourJ,
                                                                      neighbourK),
                                               *neigh_proc_id);

                  // If the neighbour is on this rank, ignore it.
                  if (topology::NetworkTopology::Instance()->GetLocalRank() == *neigh_proc_id)
                  {
                    continue;
                  }

                  if (neighbouringProcessors.count(*neigh_proc_id) == 0)
                  {
                    NeighbouringProcessor newProc(*neigh_proc_id);
                    neighbouringProcessors[*neigh_proc_id] = newProc;
                  }
                }
              }
            }
          }
          while (siteTraverser.TraverseOne());
        }
        while (blockTraverser.TraverseOne());

        // Iterate over the blocks, updating the value of the counter variable wherever
        // there is velocity field data.
        for (site_t block = 0; block < latDat.GetBlockCount(); block++)
        {
          if (velocityField[block].empty())
          {
            continue;
          }

          if (latDat.GetBlock(block)->site_data == NULL)
          {
            continue;
          }

          // Update the site id on each velocity field unit as required.
          for (site_t localSiteId = 0; localSiteId < latDat.GetSitesPerBlockVolumeUnit(); localSiteId++)
          {
            velocityField[block][localSiteId].site_id
                = latDat.GetBlock(block)->site_data[localSiteId];
          }
        }
      }

      bool VelocityField::BlockContainsData(size_t blockNumber) const
      {
        return !velocityField[blockNumber].empty();
      }

      VelocitySiteData& VelocityField::GetSiteData(site_t blockNumber, site_t siteNumber)
      {
        return velocityField[blockNumber][siteNumber];
      }

      // Returns the velocity site data for a given index, or NULL if the index isn't valid / has
      // no data.
      VelocitySiteData* VelocityField::GetVelocitySiteData(const geometry::LatticeData& latDat,
                                                           const util::Vector3D<site_t>& location)
      {
        if (!latDat.IsValidLatticeSite(location.x, location.y, location.z))
        {
          return NULL;
        }

        // TODO this stuff should be encapsulated in the LatticeData
        site_t blockI = location.x >> latDat.GetLog2BlockSize();
        site_t blockJ = location.y >> latDat.GetLog2BlockSize();
        site_t blockK = location.z >> latDat.GetLog2BlockSize();

        site_t block_id = latDat.GetBlockIdFromBlockCoords(blockI, blockJ, blockK);

        if (!BlockContainsData(static_cast<size_t> (block_id)))
        {
          return NULL;
        }

        site_t localSiteI = location.x - (blockI << latDat.GetLog2BlockSize());
        site_t localSiteJ = location.y - (blockJ << latDat.GetLog2BlockSize());
        site_t localSiteK = location.z - (blockK << latDat.GetLog2BlockSize());

        site_t site_id = ( ( (localSiteI << latDat.GetLog2BlockSize()) + localSiteJ)
            << latDat.GetLog2BlockSize()) + localSiteK;

        return &GetSiteData(block_id, site_id);
      }

      // Function to initialise the velocity field at given coordinates.
      void VelocityField::InitializeVelocityFieldBlock(const geometry::LatticeData& latDat,
                                                       const util::Vector3D<site_t> location,
                                                       const proc_t proc_id)
      {
        // TODO this stuff should be encapsulated in the LatticeData
        site_t blockI = location.x >> latDat.GetLog2BlockSize();
        site_t blockJ = location.y >> latDat.GetLog2BlockSize();
        site_t blockK = location.z >> latDat.GetLog2BlockSize();

        site_t blockId = latDat.GetBlockIdFromBlockCoords(blockI, blockJ, blockK);

        if (!BlockContainsData(blockId))
        {
          velocityField[blockId]
              = std::vector<VelocitySiteData>(latDat.GetSitesPerBlockVolumeUnit());
        }

        site_t localSiteI = location.x - (blockI << latDat.GetLog2BlockSize());
        site_t localSiteJ = location.y - (blockJ << latDat.GetLog2BlockSize());
        site_t localSiteK = location.z - (blockK << latDat.GetLog2BlockSize());

        site_t localSiteId = ( ( (localSiteI << latDat.GetLog2BlockSize()) + localSiteJ)
            << latDat.GetLog2BlockSize()) + localSiteK;
        velocityField[blockId][localSiteId].proc_id = proc_id;
      }

      // Populate the matrix v with all the velocity field data at each index.
      // Returns true if this the area resides entirely on this core.
      void VelocityField::GetVelocityFieldAroundPoint(const util::Vector3D<site_t> location,
                                                      const geometry::LatticeData& latDat,
                                                      util::Vector3D<float> localVelocityField[2][2][2])
      {
        for (int unitGridI = 0; unitGridI <= 1; ++unitGridI)
        {
          site_t neighbourI = location.x + unitGridI;

          for (int unitGridJ = 0; unitGridJ <= 1; ++unitGridJ)
          {
            site_t neighbourJ = location.y + unitGridJ;

            for (int unitGridK = 0; unitGridK <= 1; ++unitGridK)
            {
              site_t neighbourK = location.z + unitGridK;

              if (!latDat.IsValidLatticeSite(neighbourI, neighbourJ, neighbourK))
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                localVelocityField[unitGridI][unitGridJ][unitGridK] = util::Vector3D<float>::Zero();
                continue;
              }

              VelocitySiteData *vel_site_data_p =
                  GetVelocitySiteData(latDat,
                                      util::Vector3D<site_t>(neighbourI, neighbourJ, neighbourK));

              if (vel_site_data_p->proc_id == -1)
              {
                // it is a solid site and the velocity is
                // assumed to be zero
                localVelocityField[unitGridI][unitGridJ][unitGridK] = util::Vector3D<float>::Zero();
                continue;
              }

              if (vel_site_data_p->counter != counter)
              {
                UpdateLocalField(vel_site_data_p, latDat);
              }

              localVelocityField[unitGridI][unitGridJ][unitGridK] = vel_site_data_p->velocity;
            }
          }
        }
      }

      bool VelocityField::NeededFromNeighbour(const util::Vector3D<site_t> location,
                                              const geometry::LatticeData& latDat,
                                              proc_t* sourceProcessor)
      {
        const proc_t thisRank = topology::NetworkTopology::Instance()->GetLocalRank();

        if (!latDat.IsValidLatticeSite(location.x, location.y, location.z))
        {
          return false;
        }

        VelocitySiteData *vel_site_data_p = GetVelocitySiteData(latDat, location);

        if (vel_site_data_p->proc_id == -1 || vel_site_data_p->proc_id == thisRank
            || vel_site_data_p->counter == counter)
        {
          return false;
        }

        vel_site_data_p->counter = counter;
        *sourceProcessor = vel_site_data_p->proc_id;

        return true;
      }

      void VelocityField::UpdateLocalField(const util::Vector3D<site_t>& position,
                                           const geometry::LatticeData& latDat)
      {
        VelocitySiteData *localVelocitySiteData = GetVelocitySiteData(latDat, position);

        if (log::Logger::ShouldDisplay<log::Debug>())
        {
          if (topology::NetworkTopology::Instance()->GetLocalRank()
              != localVelocitySiteData->proc_id)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Got a request for velocity data "
              "that actually seems to be on rank %i", localVelocitySiteData->proc_id);
          }
        }

        UpdateLocalField(localVelocitySiteData, latDat);
      }

      void VelocityField::UpdateLocalField(VelocitySiteData* localVelocitySiteData,
                                           const geometry::LatticeData& latDat)
      {
        // the local counter is set equal to the global one
        // and the local velocity is calculated
        localVelocitySiteData->counter = counter;
        distribn_t density, vx, vy, vz;

        D3Q15::CalculateDensityAndVelocity(latDat.GetFOld(localVelocitySiteData->site_id
                                               * D3Q15::NUMVECTORS),
                                           density,
                                           vx,
                                           vy,
                                           vz);

        localVelocitySiteData->velocity.x = (float) (vx / density);
        localVelocitySiteData->velocity.y = (float) (vy / density);
        localVelocitySiteData->velocity.z = (float) (vz / density);
      }

      // Interpolates a velocity field to get the velocity at the position of a particle.
      util::Vector3D<float> VelocityField::InterpolateVelocityForPoint(const util::Vector3D<float> point,
                                                                       const util::Vector3D<float> localVelocityField[2][2][2]) const
      {
        float dummy;

        // Get the fractional parts of each of x, y, z
        float dx = modff(point.x, &dummy);
        float dy = modff(point.y, &dummy);
        float dz = modff(point.z, &dummy);

        util::Vector3D<float> yInterpolatedVelocity[2];

        for (int unitX = 0; unitX <= 1; unitX++)
        {
          util::Vector3D<float> zInterpolatedVelocityForThisX[2];

          for (int unitY = 0; unitY <= 1; unitY++)
          {
            zInterpolatedVelocityForThisX[unitY] = localVelocityField[unitX][unitY][0] * (1.F - dz)
                + localVelocityField[unitX][unitY][1] * dz;
          }

          yInterpolatedVelocity[unitX] = zInterpolatedVelocityForThisX[0] * (1.F - dy)
              + zInterpolatedVelocityForThisX[1] * dy;
        }

        return yInterpolatedVelocity[0] * (1.F - dx) + yInterpolatedVelocity[1] * dx;
      }

      void VelocityField::InvalidateAllCalculatedVelocities()
      {
        ++counter;
      }

    }
  }
}
