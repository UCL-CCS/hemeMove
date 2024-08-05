// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/BodyForces.h"
#include "colloids/BodyForceExamples.h"
#include "colloids/GraviticBodyForce.h"

namespace hemelb::colloids
{
    std::map<std::string, const BodyForce* const > BodyForces::bodyForces;
    std::map<site_t, LatticeForceVector> BodyForces::forceForEachSite;

    const void BodyForces::InitBodyForces(io::xml::Document& xml)
    {
      std::map<std::string, BodyForceFactory_Create> mapForceGenerators;
      mapForceGenerators["gravitic"] = & (GraviticBodyForceFactory::Create);
      mapForceGenerators["constant"] = & (ConstantBodyForceFactory::Create);
      mapForceGenerators["inv_r_sq"] = & (RadialBodyForceFactory::Create);

      io::xml::Element colloidsBodyForcesNode =
          xml.GetRoot().GetChildOrThrow("colloids").GetChildOrThrow("bodyForces");

      for (auto& [forceClass, createFunction]: mapForceGenerators)
      {
        for (auto forceNode: colloidsBodyForcesNode.Children(forceClass))
        {
          auto forceName = forceNode.GetAttributeOrThrow("forceName");
          BodyForce* nextForce = createFunction(forceNode);
          BodyForces::bodyForces.insert(std::make_pair(forceName, nextForce));
        }
      }
    }

    const LatticeForceVector BodyForces::GetBodyForcesForParticle(const Particle& particle)
    {
      auto totalForce = LatticeForceVector::Zero();
      for (auto iter = bodyForces.begin(); iter != bodyForces.end(); ++iter) {
        totalForce += iter->second->GetForceForParticle(particle);
      }
      return totalForce;
    }

}
