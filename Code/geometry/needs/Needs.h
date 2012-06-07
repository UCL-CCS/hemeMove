#ifndef HEMELB_GEOMETRY_NEEDS_NEEDS_H
#define HEMELB_GEOMETRY_NEEDS_NEEDS_H
#include <vector>
#include "net/net.h"
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace geometry
  {
    /*
     JH Note: Class created to make the new needs communication strategy in #142 testable.
     */
    /***
     *  Class defining HemeLB needs communication
     Used by geometry reader to know where to send which blocks.
     @tparam Net Class implementing the Net Communication Protocol, used to share information.
     */
    template<class Net> class NeedsBase
    {
      public:
        /***
         * Constructor for Needs manager.
         * @param BlockCount Count of blocks
         * @param readBlock Temporary input - which cores need which blocks, as an array of booleans.
         * @param areadingGroupSize Number of cores to use for reading blocks
         * @param anet Instance of Net Communication Protocol -- typically just Net(comm)
         * @param comm Communicator used for decomposition topology
         * @param rank Rank in decomposition topology
         * @param size Size of decomposiiton topology
         */
       NeedsBase(const site_t BlockCount,
                          const std::vector<bool>& readBlock,
                          const proc_t areadingGroupSize,
                          Net &anet,
                          bool shouldValidate); // Temporarily during the refactor, constructed just to abstract the block sharing bit

        /***
         * Which processors need a given block?
         * @param block Block number to query
         * @return Vector of ranks in the decomposition topology which need this block
         */
        const std::vector<proc_t> & ProcessorsNeedingBlock(const site_t &block) const
        {
          return procsWantingBlocksBuffer[block];
        }

        /***
         * Which core should be responsible for reading a given block? This core does not necessarily
         * require information about the block
         *
         * @param blockNumber Block number to query
         * @return Rank in the decomposition topology, for core which should read the block.
         */
        proc_t GetReadingCoreForBlock(const site_t blockNumber) const;
      private:
        std::vector<std::vector<proc_t> > procsWantingBlocksBuffer;
        Net &net;
        const topology::Communicator & communicator;
        const proc_t readingGroupSize;
        bool shouldValidate;
        void Validate(const site_t blockCount, const std::vector<bool>& readBlock);
    };
    typedef NeedsBase<net::Net> Needs;
  }
}
#endif // HEMELB_GEOMETRY_NEEDS_NEEDS_H
