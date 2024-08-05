// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GMYREADRESULT_H
#define HEMELB_GEOMETRY_GMYREADRESULT_H

#include <memory>
#include <map>
#include <vector>

#include "units.h"
#include "constants.h"
#include "geometry/GeometryBlock.h"
#include "util/numerical.h"
#include "util/Vector3D.h"

namespace hemelb::geometry
{
    namespace octree { class DistributedStore; }

    // A site description is the block's OCT index (hence U64 to
    // match the tree's choice) and the index of the site within
    // the block (could be U16, but use 64 as padding will force
    // that anyway).
    using SiteDesc = std::array<U64,2>;
    using SiteVec = std::vector<SiteDesc>;
    using MovesMap = std::map<int, SiteVec>;

    /***
     * Model of the information in a geometry file
     */
    class GmyReadResult
    {
    public:
        /**
         * Default constructor initialises internal variables
         */
        GmyReadResult(const Vec16& dimensionsInBlocks, U16 blockSize);

        //! Destructor has to be defined in a TU where LookupTree is defined
        ~GmyReadResult();
        GmyReadResult(GmyReadResult&&) = default;
        GmyReadResult& operator=(GmyReadResult&&) = default;

        /***
         * Returns the total count of blocks in the bounding box of the geometry.
         * @return count of blocks in the geometry.
         */
        inline site_t GetBlockCount() const
        {
          return blockCount;
        }

        /**
         * Returns the number of sites in each cubic block.
         * @return Number of sites in each block of the geometry.
         */
        inline site_t GetSitesPerBlock() const
        {
          return sitesPerBlock;
        }

        /***
         * Get the i.d. of a block, i.e. the one-d coordinate, from the three-d coordinate.
         * @todo Use this to replace domain_type::GetBlockGmyIdxFromBlockCoords
         */
        inline site_t GetBlockIdFromBlockCoordinates(U16 blockI, U16 blockJ, U16 blockK) const
        {
          return (blockI * dimensionsInBlocks.y() + blockJ) * dimensionsInBlocks.z() + blockK;
        }

        inline U64 GetBlockIdFromBlockCoordinates(Vec16 ijk) const {
            auto IJK = ijk.as<U64>();
            return (IJK[0] * dimensionsInBlocks[1] + IJK[1]) * dimensionsInBlocks[2] + IJK[2];
        }

        /***
         * Get the coordinates of a block from a block id.
         */
        inline Vec16 GetBlockCoordinatesFromBlockId(site_t blockId) const
        {
          site_t blockZ = blockId % dimensionsInBlocks.z();
          site_t remainder = blockId / dimensionsInBlocks.z();
          site_t blockY = remainder % dimensionsInBlocks.y();
          site_t blockX = remainder / dimensionsInBlocks.y();
          return Vec16(blockX, blockY, blockZ);
        }

        /***
         * Get the i.d. of a site, i.e. the one-d coordinate, from the three-d coordinate.
         */
        inline site_t GetSiteIdFromSiteCoordinates(site_t siteI, site_t siteJ, site_t siteK) const
        {
          return (siteI * blockSize + siteJ) * blockSize + siteK;
        }

        /**
         * True if the given block coordinates are within the geometry bounding-box.
         */
        inline bool AreBlockCoordinatesValid(const Vec16& blockCoords) const
        {
          return blockCoords.x() >= 0 && blockCoords.y() >= 0 && blockCoords.z() >= 0
              && blockCoords.x() < dimensionsInBlocks.x() && blockCoords.y() < dimensionsInBlocks.y()
              && blockCoords.z() < dimensionsInBlocks.z();
        }

        /**
         * True if the given site coordinates are within a bounding-box for a block.
         */
        inline bool AreLocalSiteCoordinatesValid(const util::Vector3D<site_t>& siteCoords) const
        {
          return siteCoords.x() >= 0 && siteCoords.y() >= 0 && siteCoords.z() >= 0
              && siteCoords.x() < blockSize && siteCoords.y() < blockSize && siteCoords.z() < blockSize;
        }

        /**
         * Get the dimensions of the bounding box in terms of blocks.
         * @return Dimensions of the bounding box in blocks.
         */
        inline Vec16 const& GetBlockDimensions() const
        {
          return dimensionsInBlocks;
        }

        /**
         * Get the number of sites along one edge of a block.
         * @return Number of sites that make up one block-length.
         */
        inline U16 GetBlockSize() const
        {
          return blockSize;
        }

        /* Find a site index, taking into account ALL lattice sites. */
        site_t FindSiteIndexInBlock(site_t fluidSiteBlock, site_t fluidSitesToPass) const;

        /* Find a site index, taking into account ONLY fluid sites. */
        site_t FindFluidSiteIndexInBlock(site_t fluidSiteBlock, site_t neighbourSiteId) const;

      private:
        Vec16 dimensionsInBlocks; //! The count of blocks in each direction
        U16 blockSize; //! Size of a block, in sites.

        site_t blockCount;
        site_t sitesPerBlock;

      public:
        std::vector<BlockReadResult> Blocks; //! Array of Block models
        std::unique_ptr<octree::DistributedStore> block_store;
    };

    // This iterator will move through the sites that make up a
    // block. On dereference it returns a pair of the current 3D
    // position and the 1D index in geometry file order.
    struct BlockSiteIterator {
        GmyReadResult const* gmy;
        Vec16 pos;
        site_t idx;

        using value_type = std::pair<Vec16, site_t>;
        inline value_type operator*() const {
            return {pos, idx};
        }

        inline BlockSiteIterator& operator++() {
            auto B = gmy->GetBlockSize();
            idx++;
            if (++pos[2] == B) {
                pos[2] = 0;
                if (++pos[1] == B) {
                    pos[1] = 0;
                    ++pos[0];
                }
            }
            return *this;
        }

        inline friend bool operator==(BlockSiteIterator const& l, BlockSiteIterator const& r) {
          return l.idx == r.idx;
        }
    };

    // Represent the range of valid sites in a block that belongs to
    // the GmyReadResult.
    struct BlockSiteRange {
        GmyReadResult const* gmy;
        auto begin() const {
            return BlockSiteIterator{gmy, {0,0,0}, 0};
        }
        auto end() const {
            auto B = gmy->GetBlockSize();
            return BlockSiteIterator{gmy, {B, B, B}, B*B*B};
        }
    };

    // Helper for the above.
    inline auto IterSitesInBlock(GmyReadResult const& gmy) {
        return BlockSiteRange{&gmy};
    }
}
#endif // HEMELB_GEOMETRY_GMYREADRESULT_H
