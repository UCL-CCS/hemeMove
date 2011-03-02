#ifndef HEMELB_LB_GLOBALLATTICEDATA_H
#define HEMELB_LB_GLOBALLATTICEDATA_H

#include "D3Q15.h"
#include "lb/LocalLatticeData.h"
#include <cstdlib>

namespace hemelb
{
  namespace lb
  {
    // Data about an element of the domain wall
    struct WallData
    {
        // estimated wall normal (if the site is close to the wall);
        double wall_nor[3];
        // cut distances along the 14 non-zero lattice vectors;
        // each one is between 0 and 1 if the surface cuts the corresponding
        // vector or is equal to "BIG_NUMBER" otherwise
        double cut_dist[D3Q15::NUMVECTORS - 1];
    };

    // Data about each global block in the lattice,
    // site_data[] is an array containing individual lattice site data
    // within a global block.
    struct BlockData
    {
        BlockData()
        {
          ProcessorRankForEachBlockSite = NULL;
          wall_data = NULL;
          site_data = NULL;
        }

        ~BlockData()
        {
          if (ProcessorRankForEachBlockSite != NULL)
            delete[] ProcessorRankForEachBlockSite;
          if (wall_data != NULL)
            delete[] wall_data;
          if (site_data != NULL)
            delete[] site_data;
        }

        // An array of the ranks on which each lattice site within the block resides.
        int *ProcessorRankForEachBlockSite;
        // Information about wall / inlet / outlet position and orientation for
        // each site.
        WallData *wall_data;
        // The "site data" for each site.
        unsigned int *site_data;
    };

    class GlobalLatticeData
    {
        friend class BlockCounter;

      public:
        void SetBasicDetails(unsigned int iBlocksX,
                             unsigned int iBlocksY,
                             unsigned int iBlocksZ,
                             unsigned int iBlockSize)
        {
          mBlocksX = iBlocksX;
          mBlocksY = iBlocksY;
          mBlocksZ = iBlocksZ;
          mBlockSize = iBlockSize;

          mSitesX = mBlocksX * mBlockSize;
          mSitesY = mBlocksY * mBlockSize;
          mSitesZ = mBlocksZ * mBlockSize;

          SitesPerBlockVolumeUnit = mBlockSize * mBlockSize * mBlockSize;

          // A shift value we'll need later = log_2(block_size)
          unsigned int i = mBlockSize;
          Log2BlockSize = 0;
          while (i > 1)
          {
            i >>= 1;
            ++Log2BlockSize;
          }

          mBlockCount = mBlocksX * mBlocksY * mBlocksZ;

          Blocks = new BlockData[mBlockCount];
        }

        unsigned int GetXSiteCount() const
        {
          return mSitesX;
        }

        unsigned int GetYSiteCount() const
        {
          return mSitesY;
        }

        unsigned int GetZSiteCount() const
        {
          return mSitesZ;
        }

        unsigned int GetXBlockCount() const
        {
          return mBlocksX;
        }

        unsigned int GetYBlockCount() const
        {
          return mBlocksY;
        }

        unsigned int GetZBlockCount() const
        {
          return mBlocksZ;
        }

        unsigned int GetBlockSize() const
        {
          return mBlockSize;
        }

        unsigned int GetBlockCount() const
        {
          return mBlockCount;
        }

        BlockData * Blocks;

        ~GlobalLatticeData()
        {
          delete[] Blocks;
        }

        // Returns the type of collision/streaming update for the fluid site
        // with data "site_data".
        unsigned int GetCollisionType(unsigned int site_data) const
        {
          unsigned int boundary_type;

          if (site_data == hemelb::lb::FLUID_TYPE)
          {
            return FLUID;
          }
          boundary_type = site_data & SITE_TYPE_MASK;

          if (boundary_type == hemelb::lb::FLUID_TYPE)
          {
            return EDGE;
          }
          if (! (site_data & PRESSURE_EDGE_MASK))
          {
            if (boundary_type == hemelb::lb::INLET_TYPE)
            {
              return INLET;
            }
            else
            {
              return OUTLET;
            }
          }
          else
          {
            if (boundary_type == hemelb::lb::INLET_TYPE)
            {
              return INLET | EDGE;
            }
            else
            {
              return OUTLET | EDGE;
            }
          }
        }

        // Function that finds the pointer to the rank on which a particular site
        // resides. If the site is in an empty block, return NULL.
        int * GetProcIdFromGlobalCoords(unsigned int iSiteI,
                                        unsigned int iSiteJ,
                                        unsigned int iSiteK) const
        {
          // If the given site location is outside the bounding box return a NULL
          // pointer.
          if (iSiteI >= mSitesX || iSiteJ >= mSitesY || iSiteK >= mSitesZ)
          {
            return NULL;
          }

          // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
          unsigned int i = iSiteI >> Log2BlockSize;
          unsigned int j = iSiteJ >> Log2BlockSize;
          unsigned int k = iSiteK >> Log2BlockSize;

          // Get the block from the block identifiers.
          BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(i, j, k)];

          // If an empty (solid) block is addressed, return a NULL pointer.
          if (lBlock->ProcessorRankForEachBlockSite == NULL)
          {
            return NULL;
          }
          else
          {
            // Find site coordinates within the block
            unsigned int ii = iSiteI - (i << Log2BlockSize);
            unsigned int jj = iSiteJ - (j << Log2BlockSize);
            unsigned int kk = iSiteK - (k << Log2BlockSize);

            // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
            // mProcessorsForEachBlock)
            return &lBlock->ProcessorRankForEachBlockSite[ ( ( (ii << Log2BlockSize) + jj)
                << Log2BlockSize) + kk];
          }
        }

        // Function that gets the index of a block from its coordinates.
        unsigned int GetBlockIdFromBlockCoords(unsigned int blockI,
                                               unsigned int blockJ,
                                               unsigned int blockK) const
        {
          // Get the block from the block identifiers.
          return (blockI * mBlocksY + blockJ) * mBlocksZ + blockK;
        }

        // Function to get a pointer to the site_data for a site.
        // If the site is in an empty block, return NULL.

        const unsigned int * GetSiteData(unsigned int iSiteI,
                                         unsigned int iSiteJ,
                                         unsigned int iSiteK) const
        {
          // If site is out of the bounding box, return NULL.
          if (iSiteI >= mSitesX || iSiteJ >= mSitesY || iSiteK >= mSitesZ)
          {
            return NULL;
          }

          // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
          unsigned int i = iSiteI >> Log2BlockSize;
          unsigned int j = iSiteJ >> Log2BlockSize;
          unsigned int k = iSiteK >> Log2BlockSize;

          // Pointer to the block
          BlockData * lBlock = &Blocks[GetBlockIdFromBlockCoords(i, j, k)];

          // if an empty (solid) block is addressed
          if (lBlock->site_data == NULL)
          {
            return NULL;
          }
          else
          {
            // Find site coordinates within the block
            unsigned int ii = iSiteI - (i << Log2BlockSize);
            unsigned int jj = iSiteJ - (j << Log2BlockSize);
            unsigned int kk = iSiteK - (k << Log2BlockSize);

            // Return pointer to site_data[site]
            return &lBlock->site_data[ ( ( (ii << Log2BlockSize) + jj) << Log2BlockSize) + kk];
          }
        }

      public:
        // TODO public temporarily, until all usages are internal to the class.
        unsigned int Log2BlockSize;
        unsigned int SitesPerBlockVolumeUnit;

      private:
        unsigned int mBlockCount;
        unsigned int mSitesX, mSitesY, mSitesZ;
        unsigned int mBlocksX, mBlocksY, mBlocksZ;
        unsigned int mBlockSize;
    };

    class BlockCounter
    {
      public:
        BlockCounter(const GlobalLatticeData* iGlobLatDat, unsigned int iStartNumber)
        {
          mBlockNumber = iStartNumber;
          mGlobLatDat = iGlobLatDat;
        }

        void operator++()
        {
          mBlockNumber++;
        }

        void operator++(int in)
        {
          mBlockNumber++;
        }

        operator int()
        {
          return mBlockNumber;
        }

        bool operator<(unsigned int iUpperLimit) const
        {
          return mBlockNumber < iUpperLimit;
        }

        int GetICoord()
        {
          return (mBlockNumber - (mBlockNumber % (mGlobLatDat->GetYBlockCount()
              * mGlobLatDat->GetZBlockCount()))) / (mGlobLatDat->GetYBlockCount()
              * mGlobLatDat->GetZBlockCount());
        }

        int GetJCoord()
        {
          int lTemp = mBlockNumber
              % (mGlobLatDat->GetYBlockCount() * mGlobLatDat->GetZBlockCount());
          return (lTemp - (lTemp % mGlobLatDat->GetZBlockCount())) / mGlobLatDat->GetZBlockCount();
        }

        int GetKCoord()
        {
          return mBlockNumber % mGlobLatDat->GetZBlockCount();
        }

        int GetICoord(int iSiteI)
        {
          return (GetICoord() << mGlobLatDat->Log2BlockSize) + iSiteI;
        }

        int GetJCoord(int iSiteJ)
        {
          return (GetJCoord() << mGlobLatDat->Log2BlockSize) + iSiteJ;
        }

        int GetKCoord(int iSiteK)
        {
          return (GetKCoord() << mGlobLatDat->Log2BlockSize) + iSiteK;
        }

      private:
        const GlobalLatticeData* mGlobLatDat;
        unsigned int mBlockNumber;
    };
  }
}

#endif /* HEMELB_LB_GLOBALLATTICEDATA_H */
