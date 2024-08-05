// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/FaderCell.h"
#include "redblood/FlowExtension.h"

namespace hemelb
{
  namespace redblood
  {
    LatticeEnergy FaderCell::operator()() const
    {
      auto const barycentre = wrappee->GetBarycentre();
      auto const energy = wrappee->Energy();
      for (auto const& extension : *iolets)
      {
        auto const weight = linearWeight(extension, barycentre);
        if (weight > 1e-12)
        {
          return energy * weight;
        }
      }
      return energy;
    }

    LatticeEnergy FaderCell::operator()(std::vector<LatticeForceVector> &forces) const
    {
      auto const barycentre = wrappee->GetBarycentre();
      auto const energy = wrappee->Energy(forces);
      for (auto const& extension : *iolets)
      {
        auto const weight = linearWeight(extension, barycentre);
        if (weight > 1e-12)
        {
          for (auto& force : forces)
          {
            force *= weight;
          }
          return energy * weight;
        }
      }
      return energy;
    }
  }
} // hemelb::redblood
