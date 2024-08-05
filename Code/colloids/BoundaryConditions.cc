// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/BoundaryConditions.h"
#include "colloids/LubricationBC.h"
#include "colloids/DeletionBC.h"
#include "geometry/Site.h"
#include "geometry/SiteData.h"

namespace hemelb::colloids
{
    std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsWall;
    std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsInlet;
    std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsOutlet;

    const geometry::Domain* BoundaryConditions::latticeData;

    const void BoundaryConditions::InitBoundaryConditions(
            const geometry::Domain* const latticeData, io::xml::Document& xml)
    {
      BoundaryConditions::latticeData = latticeData;

      std::map<std::string, BoundaryConditionFactory_Create> mapBCGenerators;
      mapBCGenerators["lubricationBC"] = & (LubricationBoundaryConditionFactory::Create);
      mapBCGenerators["deletionBC"] = & (DeletionBoundaryConditionFactory::Create);

      io::xml::Element colloidsBC =
          xml.GetRoot().GetChildOrThrow("colloids").GetChildOrThrow("boundaryConditions");

      for (auto const& [boundaryConditionClass, createFunction]: mapBCGenerators)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("*** In BoundaryConditions::InitBoundaryConditions - looking for %s BC in XML\n",
                                                      boundaryConditionClass.c_str());
        // There must be at least one BC element for each type - count here
        unsigned num_bc_type = 0;
        for (auto bcNode: colloidsBC.Children(boundaryConditionClass)) {
            ++num_bc_type;
            auto appliesTo = bcNode.GetAttributeOrThrow("appliesTo");
            BoundaryCondition* nextBC = createFunction(bcNode);
            if (appliesTo == "wall")
                BoundaryConditions::boundaryConditionsWall.push_back(nextBC);
            else if (appliesTo == "inlet")
                BoundaryConditions::boundaryConditionsInlet.push_back(nextBC);
            else if (appliesTo == "outlet")
                BoundaryConditions::boundaryConditionsOutlet.push_back(nextBC);
        }
        if (num_bc_type == 0)
            throw (Exception() << "Missing colloids boundary condition element for class "
                               << boundaryConditionClass);
      }
    }

    const bool BoundaryConditions::DoSomeThingsToParticle(const LatticeTimeStep currentTimestep,
                                                          Particle& particle)
    {
      bool keep = true;

      particle.SetLubricationVelocityAdjustment(LatticeVelocity());

      // detect collision(s)
      LatticePosition const &pos = particle.GetGlobalPosition();
      auto const siteGlobalPosition = (pos + util::Vector3D(0.5)).as<site_t>();
      proc_t procId;
      site_t localContiguousId;
      const bool isLocalFluid = latticeData->GetContiguousSiteId(siteGlobalPosition,
                                                                 procId,
                                                                 localContiguousId);
      if (pos.y() < 1.5 && pos.y() >= 0.5)
        log::Logger::Log<log::Trace, log::OnePerCore>("*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu, p.pos: {%g,%g,%g}, p.vel: {%g,%g,%g}, isLocalFluid: %s, procId: %u, localContiguousId: %lu, siteCoords: {%lu,%lu,%lu}, ownerRank: %u\n",
                                                      particle.GetParticleId(),
                                                      pos.x(),
                                                      pos.y(),
                                                      pos.z(),
                                                      particle.GetVelocity().x(),
                                                      particle.GetVelocity().y(),
                                                      particle.GetVelocity().z(),
                                                      isLocalFluid ?
                                                        "TRUE" :
                                                        "FALSE",
                                                      procId,
                                                      localContiguousId,
                                                      siteGlobalPosition.x(),
                                                      siteGlobalPosition.y(),
                                                      siteGlobalPosition.z(),
                                                      particle.GetOwnerRank());

      if (!isLocalFluid)
      {
        return keep;
      }

      auto& latticeInfo = BoundaryConditions::latticeData->GetLatticeInfo();
      const geometry::Site<const geometry::Domain> site =
          latticeData->GetSite(localContiguousId);
      const geometry::SiteData siteData = site.GetSiteData();
      const geometry::SiteType siteType = siteData.GetSiteType();
      const distribn_t* siteWallDistances = site.GetWallDistances();

      const bool isNearWall = siteData.IsWall();
      const bool isNearInlet = (siteType == geometry::INLET_TYPE);
      const bool isNearOutlet = (siteType == geometry::OUTLET_TYPE);

      // if the particle is not near a boundary then simply keep it
      if (!isNearWall && !isNearInlet && !isNearOutlet)
      {
        particle.SetDeletionMarker();
        return keep;
      }
      ////else
      log::Logger::Log<log::Trace, log::OnePerCore>("*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu, isNearWall: %s, isNearInlet: %s, isNearOutlet: %s ***\n",
                                                    particle.GetParticleId(),
                                                    isNearWall ?
                                                      "TRUE" :
                                                      "FALSE",
                                                    isNearInlet ?
                                                      "TRUE" :
                                                      "FALSE",
                                                    isNearOutlet ?
                                                      "TRUE" :
                                                      "FALSE");

      // only use lattice vectors 1 to 6 (the face-of-a-cube vectors)
      std::vector<LatticePosition> particleToWallVectors;
      for (Direction direction = 1; direction <= 6; ++direction)
      {
        // in general, this "distance" is a fraction of a non-unit lattice vector
        // however, we treat this fractional magnitude as a real lattice distance
        // because all of the face-of-a-cube lattice vectors will be unit vectors
        double thisDistance = siteWallDistances[direction - 1];

        // a negative distance to the wall from a site in any direction means that
        // the wall is further away than the nearest fluid site in that direction
        if (thisDistance < 0.0)
          continue;

        // the particle cannot be allowed to go past halfway between this site and
        // the next lattice site in this direction, because the next site is solid
        // the wall is assumed to be no further away than half the distance to the
        // solid site so that the particle never becomes nearest to a solid site
        if (thisDistance > 0.5)
          thisDistance = 0.5;

        // conversion from LatticeCoordinate to LatticePosition is done
        // auto-magically by the multiplication & its arithmetic traits
        const LatticePosition siteToWall = latticeInfo.GetVector(direction) * thisDistance;
        const LatticePosition particleToSite = (LatticePosition) siteGlobalPosition
                                               - pos;

        // particleToWall = siteToWall + projection of particleToSite in the siteToWall direction
        const LatticePosition particleToWallVector = siteToWall
            + siteToWall.GetNormalised() * Dot(siteToWall.GetNormalised(), particleToSite);

        log::Logger::Log<log::Trace, log::OnePerCore>("*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu, siteToWall: {%g,%g,%g}, particleToSite: {%g,%g,%g}, particleToWall: {%g,%g,%g}\n",
                                                      particle.GetParticleId(),
                                                      siteToWall.x(),
                                                      siteToWall.y(),
                                                      siteToWall.z(),
                                                      particleToSite.x(),
                                                      particleToSite.y(),
                                                      particleToSite.z(),
                                                      particleToWallVector.x(),
                                                      particleToWallVector.y(),
                                                      particleToWallVector.z());

        particleToWallVectors.push_back(particleToWallVector);
      }

      if (isNearWall)
        for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsWall.begin();
            iter != boundaryConditionsWall.end(); iter++)
        {
          BoundaryCondition& boundaryCondition = ** (iter);
          keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
        }

      if (isNearInlet)
        for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsInlet.begin();
            iter != boundaryConditionsInlet.end(); iter++)
        {
          BoundaryCondition& boundaryCondition = ** (iter);
          keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
        }

      if (isNearOutlet)
        for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsOutlet.begin();
            iter != boundaryConditionsOutlet.end(); iter++)
        {
          BoundaryCondition& boundaryCondition = ** (iter);
          keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
        }

      if (keep)
        particle.SetDeletionMarker();
      else
      {
        particle.SetDeletionMarker(currentTimestep);
        log::Logger::Log<log::Trace, log::OnePerCore>("*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu - attempting to set markedForDeletion to %lu (value actually becomes: %lu)\n",
                                                      particle.GetParticleId(),
                                                      currentTimestep,
                                                      particle.GetDeletionMarker());
      }

      return keep;
    }

}
