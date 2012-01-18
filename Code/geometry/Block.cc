#include "constants.h"
#include "geometry/Block.h"

namespace hemelb
{
  namespace geometry
  {
    Block::Block()
    {
    }

    Block::Block(site_t sitesPerBlock) :
        processorRankForEachBlockSite(sitesPerBlock, BIG_NUMBER2), localContiguousIndex(sitesPerBlock,
                                                                                        BIG_NUMBER3)
    {
    }

    Block::~Block()
    {
    }

    bool Block::IsEmpty() const
    {
      return localContiguousIndex.empty();
    }

    proc_t Block::GetProcessorRankForSite(site_t localSiteIndex) const
    {
      return processorRankForEachBlockSite[localSiteIndex];
    }

    site_t Block::GetLocalContiguousIndexForSite(site_t localSiteIndex) const
    {
      return localContiguousIndex[localSiteIndex];
    }

    void Block::SetProcessorRankForSite(site_t localSiteIndex, proc_t rank)
    {
      processorRankForEachBlockSite[localSiteIndex] = rank;
    }

    void Block::SetLocalContiguousIndexForSite(site_t localSiteIndex, site_t contiguousIndex)
    {
      localContiguousIndex[localSiteIndex] = contiguousIndex;
    }

  }
}
