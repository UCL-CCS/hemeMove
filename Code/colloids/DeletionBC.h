// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_DELETIONBC_H
#define HEMELB_COLLOIDS_DELETIONBC_H

#include "colloids/BoundaryConditions.h"
#include "configuration/SimConfig.h"
#include "configuration/SimConfigReader.h"

namespace hemelb
{
  namespace colloids
  {
    class DeletionBC : public BoundaryCondition
    {
      public:
        static BoundaryCondition* ReadFromXml(io::xml::Element& xml)
        {
          // activationDistance defaults to 1/4 lattice unit *beyond* the boundary
          LatticeDistance activationDistance = -0.25;
          // TODO: Convert to lattice units
          io::xml::Element actDistEl = xml.GetChildOrNull("activationDistance");
          if (actDistEl != io::xml::Element::Missing())
            configuration::GetDimensionalValue(actDistEl, "lattice", activationDistance);

          return new DeletionBC(activationDistance);
        }

        virtual const bool DoSomethingToParticle(
            Particle& particle, const std::vector<LatticePosition> particleToWallVectors)
        {
          // TODO: does not do *beyond* just *within* activation distance of boundary
          //LatticeDistance distance = wallNormal.GetMagnitudeSquared();
          log::Logger::Log<log::Trace, log::OnePerCore>("*** In DeletionBC::DoSomethingToParticle for particleId: %lu ***\n",
                                                        particle.GetParticleId());
          return false; //distance < (activationDistance * activationDistance);
        }

        virtual const std::vector<Particle> CreateNewParticles()
        {
          return std::vector<Particle>();
        }

      protected:
        DeletionBC(LatticeDistance activationDistance) :
            activationDistance(activationDistance)
        {
        }
        ;

        LatticeDistance activationDistance;
    };

    class DeletionBoundaryConditionFactory : public BoundaryConditionFactory<DeletionBC>
    {
    };
  }
}
#endif /* HEMELB_COLLOIDS_DELETIONBC_H */
