// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLCONTROLLER_H
#define HEMELB_REDBLOOD_CELLCONTROLLER_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "redblood/CellArmy.h"
#include "log/Logger.h"

namespace hemelb::redblood
{
    //! \brief Federates the cells together so we can apply ops simultaneously
    //! \tparam TRAITS holds type of kernel and stencil
    template<class TRAITS>
    class CellController : public CellArmy<TRAITS>,
                           public net::IteratedAction
    {
      public:
        using Traits = TRAITS;
        using CellArmy<TRAITS>::CellArmy;

        void RequestComms() override
        {
          using namespace log;
          Logger::Log<Debug, Singleton>("Cell insertion");
          CellArmy<TRAITS>::CallCellInsertion();
          if constexpr (build_info::USE_KRUEGER_ORDERING) {
              Logger::Log<Debug, Singleton>("Cell interaction with fluid");
              CellArmy<TRAITS>::Cell2FluidInteractions();
              Logger::Log<Debug, Singleton>("Fluid interaction with cells");
              CellArmy<TRAITS>::Fluid2CellInteractions();
          } else {
              Logger::Log<Debug, Singleton>("Fluid interaction with cells");
              CellArmy<TRAITS>::Fluid2CellInteractions();
              Logger::Log<Debug, Singleton>("Cell interaction with fluid");
              CellArmy<TRAITS>::Cell2FluidInteractions();
          }
        }

        void EndIteration() override
        {
          using namespace log;
          Logger::Log<Debug, Singleton>("Checking whether cells have reached outlets");
          CellArmy<TRAITS>::CellRemoval();
          Logger::Log<Debug, Singleton>("Notify cell listeners");
          CellArmy<TRAITS>::NotifyCellChangeListeners();
        }
    };
  }

#endif
