#include "colloids/ColloidController.h"
#include "geometry/BlockTraverser.h"
#include "geometry/SiteTraverser.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace colloids
  {

    // destructor
    ColloidController::~ColloidController()
    {
    }

    // constructor - called by SimulationMaster::Initialise()
    ColloidController::ColloidController(net::Net* net,
                                         geometry::LatticeData* latDatLBM,
                                         geometry::GeometryReadResult* gmyResult) :
            mNet(net), mLatDat(latDatLBM),
            mLocalRank(topology::NetworkTopology::Instance()->GetLocalRank())
    {
      // The neighbourhood used here is different to the latticeInfo used to create latDatLBM
      // The portion of the geometry input file that was read in by this proc, i.e. gmyResult
      // contains more information than was used when creating the LB lattice, i.e. latDatLBM
      // so, we traverse mLatDat to find local fluid sites but then get neighbour information
      // from the geometry file using a neighbour lattice definition appropriate for colloids

      // get the description of the colloid neighbourhood (as a vector of Vector3D of site_t)
      neighbourhood_t neighbourhood = GetNeighbourhoodVectors((site_t)2);

      // determine information about neighbour sites and processors for all local fluid sites
      InitialiseNeighbourList(gmyResult, neighbourhood);
    }

    void ColloidController::InitialiseNeighbourList(
                                         geometry::GeometryReadResult* gmyResult,
                                         neighbourhood_t neighbourhood)
    {
      // PLAN
      // foreach block in gmyResult (i.e. each block that may have been read from the input file)
      //   if block has sites (i.e. if this process _has_ read this block in from the input file)
      //     foreach site in block (i.e. each site, which may be local or remote, fluid or solid)
      //       if site is local (i.e. if the targetProcessor for this site is equal to localRank)
      //         foreach neighbour of site (i.e. follow each direction vector in the LatticeInfo)
      //           if neighbour is valid (i.e. neighbour site is within the geometry & not solid)
      //             if neighbour is remote (i.e. targetProcessor for neighbour is not localRank)
      //               if new neighbourRank (i.e. neighbourRank is not already in neighbourRanks)
      //                 add the targetProcessor of the neighbour site to our neighbourRanks list

      // foreach block in geometry
      for (geometry::BlockTraverser blockTraverser(*mLatDat);
           blockTraverser.CurrentLocationValid();
           blockTraverser.TraverseOne())
      {
        util::Vector3D<site_t> globalLocationForBlock =
              blockTraverser.GetCurrentLocation() * gmyResult->blockSize;

        // if block has sites
        site_t blockId = blockTraverser.GetCurrentIndex();
        if (gmyResult->Blocks[blockId].Sites.size() == 0)
          continue;

        // foreach site in block
        for (geometry::SiteTraverser siteTraverser = blockTraverser.GetSiteTraverser();
             siteTraverser.CurrentLocationValid();
             siteTraverser.TraverseOne())
        {
          util::Vector3D<site_t> globalLocationForSite =
                globalLocationForBlock + siteTraverser.GetCurrentLocation();

          // if site is local
          site_t siteId = siteTraverser.GetCurrentIndex();
          if (gmyResult->Blocks[blockId].Sites[siteId].targetProcessor != mLocalRank)
            continue;

          // foreach neighbour of site
          for (neighbourhood_t::iterator itDirectionVector = neighbourhood.begin();
               itDirectionVector != neighbourhood.end();
               itDirectionVector++)
          {
            util::Vector3D<site_t> globalLocationForNeighbourSite = globalLocationForSite +
                                                                    *itDirectionVector;

            // if neighbour is valid
            site_t neighbourBlockId, neighbourSiteId;
            proc_t neighbourRank;
            bool isValid = GetLocalInformationForGlobalSite(
                  gmyResult, globalLocationForNeighbourSite,
                  &neighbourBlockId, &neighbourSiteId, &neighbourRank);

            // if neighbour is remote
            if (!isValid || neighbourRank == mLocalRank)
              continue;

            // if new neighbourRank
            int addedAlready = count(mNeighbourProcessors.begin(),
                                     mNeighbourProcessors.end(),
                                     neighbourRank);
            if (addedAlready != 0)
              continue;

            // add the targetProcessor of the neighbour site to our neighbourRanks list
            mNeighbourProcessors.push_back(neighbourRank);

            // debug message to verify the neighbour list is the same as produced by LatticeData
            if (log::Logger::ShouldDisplay<log::Debug>())
              log::Logger::Log<log::Info, log::OnePerCore>(
                  "CC: added %i as neighbour for %i because site %i in block %i is neighbour to site %i in block %i in direction (%i,%i,%i)\n",
                  (int)neighbourRank, (int)mLocalRank,
                  (int)neighbourSiteId, (int)neighbourBlockId,
                  (int)siteId, (int)blockId,
                  (*itDirectionVector).x, (*itDirectionVector).y, (*itDirectionVector).z);

          } // end for itDirectionVector
        } // end for siteTraverser
      } // end for blockTraverser

    }

    //DJH// this function should probably be in geometry::ReadResult
    bool ColloidController::GetLocalInformationForGlobalSite(
                                      geometry::GeometryReadResult* gmyResult,
                                      util::Vector3D<site_t> globalLocationForSite,
                                      site_t* blockIdForSite,
                                      site_t* localSiteIdForSite,
                                      proc_t* ownerRankForSite)
    {
      // check for global location being outside the simulation entirely
      if (globalLocationForSite.x < (site_t)0 ||
          globalLocationForSite.y < (site_t)0 ||
          globalLocationForSite.z < (site_t)0 ||
          globalLocationForSite.x >= gmyResult->blocks.x * gmyResult->blockSize ||
          globalLocationForSite.y >= gmyResult->blocks.y * gmyResult->blockSize ||
          globalLocationForSite.z >= gmyResult->blocks.z * gmyResult->blockSize )
        return false;

      // obtain block information (3D location vector and 1D id number) for the site
      util::Vector3D<site_t> blockLocationForSite = globalLocationForSite / gmyResult->blockSize;
      *blockIdForSite = gmyResult->GetBlockIdFromBlockCoordinates(blockLocationForSite.x,
                                                                  blockLocationForSite.y,
                                                                  blockLocationForSite.z);

      // if the block does not contain any sites then return invalid
      if (gmyResult->Blocks[*blockIdForSite].Sites.size() == 0)
        return false;

      // obtain site information (3D location vector and 1D id number)
      // note: these are both local to the block that contains the site
      util::Vector3D<site_t> localSiteLocation = globalLocationForSite % gmyResult->blockSize;
      *localSiteIdForSite = gmyResult->GetSiteIdFromSiteCoordinates(localSiteLocation.x,
                                                                    localSiteLocation.y,
                                                                    localSiteLocation.z);

      // obtain the rank of the processor responsible for simulating the fluid at this site
      *ownerRankForSite = gmyResult->Blocks[*blockIdForSite].Sites[*localSiteIdForSite].targetProcessor;

      // if the rank is BIG_NUMBER2 then the site is solid not fluid so return invalid
      if (*ownerRankForSite == BIG_NUMBER2)
        return false;

      // all requested information obtained and validated so return true
      return true;
    }

    // generate a vector of Vector3D objects that describe the neighbourhood
    // produces a relative vector two all sites within distance site units in all 3 directions
    // examples: if distance==1 then the vectors will describe D3Q27 lattice pattern
    //           if distance==2 then the vectors will describe a 5x5 cube pattern
    ColloidController::neighbourhood_t ColloidController::GetNeighbourhoodVectors(site_t distance)
    {
      neighbourhood_t vectors;

      for (site_t xAdj = -distance; xAdj <= distance; xAdj++)
        for (site_t yAdj = -distance; yAdj <= distance; yAdj++)
          for (site_t zAdj = -distance; zAdj <= distance; zAdj++)
          {
            vectors.push_back(util::Vector3D<site_t>(xAdj, yAdj, zAdj));
          }

      return vectors;
    }

  }
}
