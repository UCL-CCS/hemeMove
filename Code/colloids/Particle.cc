// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/Particle.h"
#include "colloids/BodyForces.h"
#include "geometry/Domain.h"
#include "lb/LbmParameters.h"
#include "log/Logger.h"
#include "units.h"

namespace hemelb::colloids
{
    Particle::Particle(const geometry::Domain& latDatLBM,
                       const hemelb::lb::LbmParameters *lbmParams, io::xml::Element& xml) :
        PersistedParticle(xml), lbmParams(lbmParams)
    {
      // updating position with zero velocity and zero body force is necessary
      // because of the side-effect that sets owner rank from the new position
      ownerRank = SITE_OR_BLOCK_SOLID;
      velocity *= 0.0;
      bodyForces *= 0.0;
      UpdatePosition(latDatLBM);

      OutputInformation();
      if (log::Logger::ShouldDisplay<log::Trace>())
        log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::ctor, id: %i, a0: %g, ah: %g, position: {%g,%g,%g}\n",
                                                      particleId,
                                                      smallRadius_a0,
                                                      largeRadius_ah,
                                                      globalPosition.x(),
                                                      globalPosition.y(),
                                                      globalPosition.z());
    }

//    const bool Particle::operator<(const Particle& other) const
//    {
//      // ORDER BY isLocal, ownerRank, particleId
//      if (ownerRank == other.ownerRank)
//        return (particleId < other.particleId);
//      else if (ownerRank == net::IOCommunicator::Instance()->Rank())
//        return true;
//      else if (other.ownerRank == net::IOCommunicator::Instance()->Rank())
//        return false;
//      else
//        return (ownerRank < other.ownerRank);
//    }

    bool Particle::IsOwnerRankKnown(
        std::map<proc_t, std::pair<unsigned int, unsigned int> > const& map) const
    {
      return map.count(ownerRank) > 0;
    }

    bool Particle::IsReadyToBeDeleted() const
    {
      return markedForDeletionTimestep < lastCheckpointTimestep;
    }

    void Particle::OutputInformation() const
    {
      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::OutputInformation, id: %i, owner: %i, drag %g, mass %g, position: {%g,%g,%g}, velocity: {%g,%g,%g}, bodyForces: {%g,%g,%g}\n",
                                                    particleId,
                                                    ownerRank,
                                                    CalculateDragCoefficient(),
                                                    mass,
                                                    globalPosition.x(),
                                                    globalPosition.y(),
                                                    globalPosition.z(),
                                                    velocity.x(),
                                                    velocity.y(),
                                                    velocity.z(),
                                                    bodyForces.x(),
                                                    bodyForces.y(),
                                                    bodyForces.z());
    }

    // this is the exact size that the xdr data produced for this particle will occupy
    // 10 fields * 8 bytes-per-field = 80 bytes, if velocity is included in the output
    // 7 fields * 8 bytes-per-field = 56 bytes, when transient fields are not included

    void Particle::WriteToStream(const LatticeTimeStep currentTimestep,
                                       io::Writer& writer)
    {
      lastCheckpointTimestep = currentTimestep;

      writer << (uint64_t) ownerRank;
      writer << (uint64_t) particleId;
      writer << smallRadius_a0 << largeRadius_ah;
      writer << globalPosition;

      // if the following code line is ever uncommented
      // change io::formats::colloids::RecordLength to 80
      //writer << velocity.x << velocity.y << velocity.z;
    }

    void Particle::UpdatePosition(const geometry::Domain& latDatLBM)
    {
      // first, update the position: newPosition = oldPosition + velocity + bodyForces * drag
      // then,  update the owner rank for the particle based on its new position

      if (log::Logger::ShouldDisplay<log::Trace>())
        log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::UpdatePosition, id: %i,\nposition: {%g,%g,%g}\nvelocity: {%g,%g,%g}\nbodyForces: {%g,%g,%g}\n",
                                                      particleId,
                                                      globalPosition.x(),
                                                      globalPosition.y(),
                                                      globalPosition.z(),
                                                      velocity.x(),
                                                      velocity.y(),
                                                      velocity.z(),
                                                      bodyForces.x(),
                                                      bodyForces.y(),
                                                      bodyForces.z());

      globalPosition += GetVelocity();

      // round the global position of the particle to the nearest site coordinates
      auto const siteGlobalPosition = (globalPosition + util::Vector3D{0.5, 0.5, 0.5}).as<site_t>();
//      const util::Vector3D<site_t> siteGlobalPosition((site_t) (0.5 + globalPosition.x),
//                                                      (site_t) (0.5 + globalPosition.y),
//                                                      (site_t) (0.5 + globalPosition.z));

      // convert the site coordinates into a local site index and find owner rank
      proc_t procId = latDatLBM.GetProcIdFromGlobalCoords(siteGlobalPosition);
      isValid = (procId != SITE_OR_BLOCK_SOLID);
      if (isValid && (ownerRank != procId))
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Changing owner of particle %i from %i to %i - %s\n",
                                                      particleId,
                                                      ownerRank,
                                                      procId,
                                                      isValid ?
                                                        "valid" :
                                                        "INVALID");
        ownerRank = procId;
      }

      if (log::Logger::ShouldDisplay<log::Trace>())
        log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::UpdatePosition, id: %i, position is now: {%g,%g,%g}\n",
                                                      particleId,
                                                      globalPosition.x(),
                                                      globalPosition.y(),
                                                      globalPosition.z());
    }

    Dimensionless Particle::GetViscosity() const
    {
      // get fluid viscosity
      //return BLOOD_VISCOSITY_Pa_s;
      return (lbmParams->GetTau() - 0.5) * Cs2 * 1.0/*timestep*/* 1.0/*density*/;
    }

    Dimensionless Particle::CalculateDragCoefficient() const
    {
      // calculate the drag due to the fluid for this particle
      return (largeRadius_ah - smallRadius_a0)
          / (6.0 * PI * GetViscosity() * smallRadius_a0 * largeRadius_ah);
    }

    void Particle::CalculateBodyForces()
    {
      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateBodyForces, id: %i, position: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x(),
                                                    globalPosition.y(),
                                                    globalPosition.z());

      // delegate the calculation of body forces to the BodyForces class
      bodyForces = BodyForces::GetBodyForcesForParticle(*this);

      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateBodyForces, id: %i, position: {%g,%g,%g}, bodyForces: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x(),
                                                    globalPosition.y(),
                                                    globalPosition.z(),
                                                    bodyForces.x(),
                                                    bodyForces.y(),
                                                    bodyForces.z());
    }

    /** modified dirac delta function according to Peskin */
    Dimensionless diracOperation(const LatticePosition& relativePosition)
    {
      Dimensionless delta = 1.0;
      for (int xyz = 0; xyz < 3; xyz++)
      {
        const LatticeDistance rmod = fabs(relativePosition[xyz]);

        if (rmod <= 1.0)
          delta *= 0.125 * (3.0 - 2.0 * rmod + sqrt(1.0 + 4.0 * rmod - 4.0 * rmod * rmod));
        else if (rmod <= 2.0)
          delta *= 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
        else
          delta = 0.0;

      }
      return delta;
    }

    void Particle::CalculateFeedbackForces(const geometry::Domain& latDatLBM) const
    //lb::MacroscopicPropertyCache& propertyCache) const
    {
      /** CalculateFeedbackForces
       *    For each local neighbour lattice site
       *    - calculate the feedback force on each neighbour lattice site
       *    - set feedback force values into the body forces cache object
       */

      log::Logger::Log<log::Debug, log::OnePerCore>("In colloids::Particle::CalculateFeedbackForces, id: %i, position: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x(),
                                                    globalPosition.y(),
                                                    globalPosition.z());

      // determine the global coordinates of the next neighbour site:
      // nested loop - x, y, z directions semi-open interval [-2, +2)
      for (site_t x = site_t(globalPosition.x()) - 1; x < site_t(globalPosition.x()) + 3; x++)
        for (site_t y = site_t(globalPosition.y()) - 1; y < site_t(globalPosition.y()) + 3; y++)
          for (site_t z = site_t(globalPosition.z()) - 1; z < site_t(globalPosition.z()) + 3; z++)
          {
            const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

            log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateFeedbackForces, particleId: %i, position: {%g,%g,%g}, siteGlobalPosition: {%i,%i,%i}, getting siteId ...\n",
                                                          particleId,
                                                          globalPosition.x(),
                                                          globalPosition.y(),
                                                          globalPosition.z(),
                                                          siteGlobalPosition.x(),
                                                          siteGlobalPosition.y(),
                                                          siteGlobalPosition.z());

            // convert the global coordinates of the site into a local site index
            proc_t procId;
            site_t siteId;
            bool isSiteValid = latDatLBM.GetContiguousSiteId(siteGlobalPosition, procId, siteId);
            bool isSiteLocal = (procId == latDatLBM.GetLocalRank());

            log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateFeedbackForces, particleId: %i, position: {%g,%g,%g}, siteGlobalPosition: {%i,%i,%i}, got siteId: %i, isSiteValid: %s, isSiteLocal: %s, procId: %i\n",
                                                          particleId,
                                                          globalPosition.x(),
                                                          globalPosition.y(),
                                                          globalPosition.z(),
                                                          siteGlobalPosition.x(),
                                                          siteGlobalPosition.y(),
                                                          siteGlobalPosition.z(),
                                                          siteId,
                                                          isSiteValid ?
                                                            "TRUE" :
                                                            "FALSE",
                                                          isSiteLocal ?
                                                            "TRUE" :
                                                            "FALSE",
                                                          procId);

            /** TODO: implement boundary conditions for invalid/solid sites */
            if (!isSiteValid || !isSiteLocal)
              continue;

            // calculate term of the interpolation sum
            LatticePosition relativePosition(siteGlobalPosition);
            relativePosition -= globalPosition;
            LatticeForceVector contribution = bodyForces * diracOperation(relativePosition);

            log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateFeedbackForces, particleId: %i, siteIndex: %i, bodyForces: {%g,%g,%g}, contribution: {%g,%g,%g}, reading cache ...\n",
                                                          particleId,
                                                          siteId,
                                                          bodyForces.x(),
                                                          bodyForces.y(),
                                                          bodyForces.z(),
                                                          contribution.x(),
                                                          contribution.y(),
                                                          contribution.z());

            // read value of force for site index from macroscopic cache
            //LatticeForceVector partialInterpolation = propertyCache.bodyForcesCache.Get(siteId) +
            LatticeForceVector partialInterpolation = BodyForces::GetBodyForcesForSiteId(siteId)
                + contribution;

            // read value of force for site index from macroscopic cache
            //propertyCache.bodyForcesCache.Put(siteId, partialInterpolation);
            BodyForces::SetBodyForcesForSiteId(siteId, partialInterpolation);

            log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateFeedbackForces, particleId: %i, siteIndex: %i, bodyForces: {%g,%g,%g}, contribution: {%g,%g,%g}, forceOnSiteSoFar: {%g,%g,%g}\n",
                                                          particleId,
                                                          siteId,
                                                          bodyForces.x(),
                                                          bodyForces.y(),
                                                          bodyForces.z(),
                                                          contribution.x(),
                                                          contribution.y(),
                                                          contribution.z(),
                                                          partialInterpolation.x(),
                                                          partialInterpolation.y(),
                                                          partialInterpolation.z());
          }

      log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateFeedbackForces, particleId: %i, bodyForces: {%g,%g,%g}, finished\n",
                                                    particleId,
                                                    bodyForces.x(),
                                                    bodyForces.y(),
                                                    bodyForces.z());
    }

    void Particle::InterpolateFluidVelocity(const geometry::Domain& latDatLBM,
                                                  const lb::MacroscopicPropertyCache& propertyCache)
    {
      /** InterpolateFluidVelocity
       *    For each local neighbour lattice site
       *    - get velocity for the site from the macroscopic cache object
       *    - calculate site's contribution to the velocity interpolation
       *    - increment particle velocity property with this contribution
       *    - will require communication to transmit remote contributions
       */

      log::Logger::Log<log::Debug, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, id: %i, position: {%g,%g,%g}\n",
                                                    particleId,
                                                    globalPosition.x(),
                                                    globalPosition.y(),
                                                    globalPosition.z());

      velocity *= 0.0;
      // determine the global coordinates of the next neighbour site:
      // nested loop - x, y, z directions semi-open interval [-2, +2)
      for (site_t x = ((site_t) globalPosition.x()) - 1; x < ((site_t) globalPosition.x()) + 3; x++)
        for (site_t y = ((site_t) globalPosition.y()) - 1; y < ((site_t) globalPosition.y()) + 3; y++)
          for (site_t z = ((site_t) globalPosition.z()) - 1; z < ((site_t) globalPosition.z()) + 3; z++)
          {
            const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

            // convert the global coordinates of the site into a local site index
            proc_t procId;
            site_t siteId;
            bool isSiteValid = latDatLBM.GetContiguousSiteId(siteGlobalPosition, procId, siteId);
            bool isSiteLocal = (procId == latDatLBM.GetLocalRank());

            if (log::Logger::ShouldDisplay<log::Trace>())
            {
              log::Logger::Log<log::Trace, log::OnePerCore>("WAIT A MINUTE 1 ...\n");

              int globalStatus = 0; // valid && local fluid
              int blockStatus = 0; // valid && non-empty
              int siteStatus = 0; // valid && fluid

              Vec16 blockCoords, localSiteCoords;
              latDatLBM.GetBlockAndLocalSiteCoords(siteGlobalPosition,
                                                   blockCoords,
                                                   localSiteCoords);
              log::Logger::Log<log::Trace, log::OnePerCore>("WAIT A MINUTE 2 ...\n");
              if (!latDatLBM.IsValidBlock(blockCoords))
                blockStatus = 1; // invalid - out of range
              else if (latDatLBM.GetBlock(blockCoords).IsEmpty())
                blockStatus = 2; // empty - entire block is solid

              log::Logger::Log<log::Trace, log::OnePerCore>("WAIT A MINUTE 3 ...\n");

              if (!latDatLBM.IsValidLatticeSite(siteGlobalPosition))
                siteStatus = 1; // invalid - out of range
              else if (blockStatus == 2)
                siteStatus = 2; // deemed solid because block is empty
              else if (blockStatus == 0)
              {
                if (latDatLBM.GetBlock(blockCoords).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)) == SITE_OR_BLOCK_SOLID)
                {
                  siteStatus = 3; // individual site is not simulated, i.e. must be solid
                }
                else if (latDatLBM.GetBlock(blockCoords).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords))
                    == latDatLBM.GetLocalRank())
                {
                  if (latDatLBM.GetBlock(blockCoords).SiteIsSolid(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)))
                    siteStatus = 4; // individual site is local but solid (should not happen?)
                }
              }
              log::Logger::Log<log::Trace, log::OnePerCore>("WAIT A MINUTE 4 ...\n");

              if (blockStatus == 1 || (siteStatus == 1))
                globalStatus = 1; // invalid - out of range
              else if (blockStatus == 2 || (siteStatus == 3) | (siteStatus == 4))
                globalStatus = 2; // solid
              else if (latDatLBM.GetBlock(blockCoords).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords))
                  != latDatLBM.GetLocalRank())
                globalStatus = 3; // remote fluid

              log::Logger::Log<log::Trace, log::OnePerCore>("ABOUT TO DO STUFF\n");
              log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteGlobalPosition: {%i,%i,%i} - %s \n",
                                                            particleId,
                                                            siteGlobalPosition.x(),
                                                            siteGlobalPosition.y(),
                                                            siteGlobalPosition.z(),
                                                            globalStatus == 0 ?
                                                              "valid && local fluid" :
                                                              (globalStatus == 1 ?
                                                                "invalid - out-of-range" :
                                                                (globalStatus == 2 ?
                                                                  "solid" :
                                                                  (globalStatus == 3 ?
                                                                    "remote fluid" :
                                                                    "unknown status"))));
              log::Logger::Log<log::Trace, log::OnePerCore>("ABOUT TO DO MORE STUFF\n");
              log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, blockCoords: {%i,%i,%i} - %s\n",
                                                            particleId,
                                                            blockCoords.x(),
                                                            blockCoords.y(),
                                                            blockCoords.z(),
                                                            blockStatus == 0 ?
                                                              "valid && non-empty" :
                                                              (blockStatus == 1 ?
                                                                "invalid - out-of-range" :
                                                                (blockStatus == 2 ?
                                                                  "empty - entire block is solid" :
                                                                  "unknown status")));
              log::Logger::Log<log::Trace, log::OnePerCore>("ABOUT TO DO EVEN MORE STUFF\n");
              log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, localSiteCoords: {%i,%i,%i} - %s\n",
                                                            particleId,
                                                            localSiteCoords.x(),
                                                            localSiteCoords.y(),
                                                            localSiteCoords.z(),
                                                            siteStatus == 0 ?
                                                              "valid && fluid" :
                                                              (siteStatus == 1 ?
                                                                "invalid - out of range" :
                                                                (siteStatus == 2 ?
                                                                  "deemed solid because block is empty" :
                                                                  (siteStatus == 3 ?
                                                                    "deemed solid because not simulated" :
                                                                    (siteStatus == 4 ?
                                                                      "local but solid (should not happen?)" :
                                                                      "unknown status")))));
//*/
              if (siteStatus == 0)
                log::Logger::Log<log::Trace, log::OnePerCore>("ABOUT TO DO ONE LAST THING\n");
              if (siteStatus == 0)
                log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, procIdForSite: %i\n",
                                                              particleId,
                                                              latDatLBM.GetBlock(blockCoords).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)));
              log::Logger::Log<log::Trace, log::OnePerCore>("DONE LOTS OF STUFF\n");

              if (!isSiteValid || !isSiteLocal)
                if (globalStatus == 0)
                {
                  log::Logger::Log<log::Trace, log::OnePerCore>("MAJOR PROBLEM! isValid: %s, isLocal: %s, procIdForSite: %i, localRank: %i\n",
                                                                isSiteValid ?
                                                                  "T" :
                                                                  "F",
                                                                isSiteLocal ?
                                                                  "T" :
                                                                  "F",
                                                                latDatLBM.GetBlock(blockCoords).GetProcessorRankForSite(latDatLBM.GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)),
                                                                latDatLBM.GetLocalRank());
                }
            }

            /** TODO: implement boundary conditions for invalid/solid sites */
            if (!isSiteValid || !isSiteLocal)
              continue;

            // read value of velocity for site index from macroscopic cache
            // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
            util::Vector3D<double> siteFluidVelocity = propertyCache.velocityCache.Get(siteId);

            // calculate term of the interpolation sum
            LatticePosition relativePosition(siteGlobalPosition);
            relativePosition -= globalPosition;
            // TODO: should be LatticeVelocity == Vector3D<LatticeSpeed> (fix as part of #437)
            util::Vector3D<double> partialInterpolation = siteFluidVelocity
                * diracOperation(relativePosition);

            // accumulate each term of the interpolation
            velocity += partialInterpolation;

            log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteIndex: %i, fluidVelocity: {%g,%g,%g}, partialInterpolation: {%g,%g,%g}, velocitySoFar: {%g,%g,%g}\n",
                                                          particleId,
                                                          siteId,
                                                          siteFluidVelocity.x(),
                                                          siteFluidVelocity.y(),
                                                          siteFluidVelocity.z(),
                                                          partialInterpolation.x(),
                                                          partialInterpolation.y(),
                                                          partialInterpolation.z(),
                                                          velocity.x(),
                                                          velocity.y(),
                                                          velocity.z());
          }
    }

}
