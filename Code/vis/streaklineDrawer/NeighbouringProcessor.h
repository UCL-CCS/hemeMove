#ifndef HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H
#define HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H

#include <vector>

#include "net/net.h"
#include "util/Vector3D.h"
#include "vis/streaklineDrawer/Particle.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      class NeighbouringProcessor
      {
        public:
          // Constructors
          NeighbouringProcessor();
          NeighbouringProcessor(proc_t iID);

          // Functions for communicating particles.
          void AddParticleToSend(const Particle& iParticle);
          bool ParticlesToBeRetrieved();
          const Particle& PopNextReceivedParticle();
          void ClearParticleSendingList();

          void ExchangeParticleCounts(net::Net& net);
          void ExchangeParticles(net::Net& net);

          // Functions for communicating velocity data.
          void AddSiteToRequestVelocityDataFor(site_t, site_t, site_t);
          site_t GetNumberOfSitesRequestedByNeighbour();
          const util::Vector3D<float>& GetReceivedVelocityField(const size_t receivedIndex) const;
          const util::Vector3D<site_t>& GetSiteCoordsBeingRequestedByNeighbour(const size_t receivedIndex) const;
          void SetVelocityFieldToSend(const size_t sendIndex,
                                      const util::Vector3D<float>& velocityFieldToSend);

          const util::Vector3D<site_t>& GetSendingSiteCoorinates(size_t sendIndex) const;
          size_t GetNumberOfSitesRequestedByThisCore() const;
          void ClearListOfRequestedSites();

          void ExchangeSiteIdCounts(net::Net& net);
          void ExchangeSiteIds(net::Net& net);
          void ExchangeVelocitiesForRequestedSites(net::Net& net);

        private:
          size_t numberOfParticlesToSend;
          std::vector<Particle> particlesToSend;

          size_t numberOfParticlesToReceive;
          std::vector<Particle> particlesToReceive;

          site_t numberOfSiteBeingRequestedByNeighbour;
          site_t numberOfSitesRequestedByThisCore;

          std::vector<util::Vector3D<site_t> > siteCoordsRequestedByThisCore;
          std::vector<util::Vector3D<site_t> > siteCoordsRequestedByNeighbour;

          std::vector<util::Vector3D<float> > velocityFieldDataForNeighbour;
          std::vector<util::Vector3D<float> > velocityFieldDataFromNeighbour;

          proc_t neighbourRank;
      };
    }
  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_NEIGHBOURINGPROCESSOR_H
