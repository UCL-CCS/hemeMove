//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_REDBLOOD_CELLCONTROLLER_H
#define HEMELB_REDBLOOD_CELLCONTROLLER_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "redblood/CellArmy.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief Federates the cells together so we can apply ops simultaneously
    //! \tparam TRAITS holds type of kernel and stencil
    template<class TRAITS>
    class CellController : public CellArmy<TRAITS>,
                           public net::IteratedAction
    {
      public:
        typedef TRAITS Traits;
        using CellArmy<TRAITS>::CellArmy;

        void RequestComms() override
        {
          using namespace log;
          Logger::Log<Debug, Singleton>("Cell insertion");
          CellArmy<TRAITS>::CallCellInsertion();
#ifdef HEMELB_USE_KRUEGER_ORDERING
          Logger::Log<Debug, Singleton>("Cell interaction with fluid");
          CellArmy<TRAITS>::Cell2FluidInteractions();
          Logger::Log<Debug, Singleton>("Fluid interaction with cells");
          CellArmy<TRAITS>::Fluid2CellInteractions();
#else
          Logger::Log<Debug, Singleton>("Fluid interaction with cells");
          CellArmy<TRAITS>::Fluid2CellInteractions();
          Logger::Log<Debug, Singleton>("Cell interaction with fluid");
          CellArmy<TRAITS>::Cell2FluidInteractions();
#endif
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
}

#endif
