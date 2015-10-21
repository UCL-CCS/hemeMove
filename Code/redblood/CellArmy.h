//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELL_ARMY_H
#define HEMELB_UNITTESTS_REDBLOOD_CELL_ARMY_H

#include <algorithm>
#include <vector>
#include <memory>
#include <iomanip>

#include "redblood/Cell.h"
#include "redblood/CellCell.h"
#include "redblood/WallCellPairIterator.h"
#include "redblood/GridAndCell.h"
#include "redblood/FlowExtension.h"
#include "geometry/LatticeData.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief Federates the cells together so we can apply ops simultaneously
    //! \tparam TRAITS holds type of kernel and stencil
    template<class TRAITS> class CellArmy
    {
      public:
        typedef typename TRAITS::Lattice Lattice;
        typedef typename TRAITS::Kernel Kernel;
        typedef typename TRAITS::Stencil Stencil;
        //! Type of callback for listening to changes to cells
        typedef std::function<void(const CellContainer &)> CellChangeListener;

        CellArmy(geometry::LatticeData &_latDat, CellContainer const &cells,
                 LatticeDistance boxsize = 10.0, LatticeDistance halo = 2.0,
                 Node2NodeForce const &cell2Cell = {0e0, 1e0, 2},
                 Node2NodeForce const &cell2Wall = {0e0, 1e0, 2} ) :
            latticeData(_latDat), cells(cells), cellDnC(cells, boxsize, halo),
            wallDnC(createWallNodeDnC<Lattice>(_latDat, boxsize, halo)),
            cell2Cell(cell2Cell), cell2Wall(cell2Wall)
        {
        }

        //! Performs fluid to lattice interactions
        void Fluid2CellInteractions();

        //! Performs lattice to fluid interactions
        void Cell2FluidInteractions();

        CellContainer::size_type size()
        {
          return cells.size();
        }

#   ifdef HEMELB_DOING_UNITTESTS
        //! Updates divide and conquer
        void updateDNC()
        {
          cellDnC.update();
        }
        CellContainer const & GetCells() const
        {
          return cells;
        }
        DivideConquerCells const & GetDNC() const
        {
          return cellDnC;
        }
#   endif

        //! Sets up call for cell insertion
        //! Called everytime CallCellInsertion is called
        void SetCellInsertion(std::function<void(CellInserter const&)> const & f)
        {
          cellInsertionCallBack = f;
        }

        //! Calls cell insertion
        void CallCellInsertion()
        {
          if (cellInsertionCallBack)
          {
            auto callback = [this](CellContainer::value_type cell)
            {
              this->AddCell(cell);
            };
            cellInsertionCallBack(callback);
          }
        }

        //! Adds a cell change listener to be notified when cell positions change
        void AddCellChangeListener(CellChangeListener const & ccl)
        {
          cellChangeListeners.push_back(ccl);
        }

        //! Invokes the callback function to output cell positions
        void NotifyCellChangeListeners()
        {
          for (CellChangeListener ccl : cellChangeListeners)
          {
            ccl(cells);
          }
        }

        //! Sets outlets within which cells disapear
        void SetOutlets(std::vector<FlowExtension> const & olets)
        {
          outlets = olets;
        }

        //! Remove cells if they have reached outlets
        void CellRemoval();

        //! Adds input cell to simulation
        void AddCell(CellContainer::value_type cell)
        {
          auto const barycenter = cell->GetBarycenter();
          log::Logger::Log<log::Debug, log::OnePerCore>(
              "Adding cell at (%f, %f, %f)",
              barycenter.x, barycenter.y, barycenter.z
          );
          cellDnC.insert(cell);
          cells.insert(cell);
        }

        //! \brief Sets cell to cell interaction forces
        //! \details Forwards arguments to Node2NodeForce constructor.
        template<class ... ARGS> void SetCell2Cell(ARGS ... args)
        {
          cell2Cell = Node2NodeForce(std::forward<ARGS>(args)...);
        }
        //! \brief Sets cell to cell interaction forces
        //! \details Forwards arguments to Node2NodeForce constructor.
        template<class ... ARGS> void SetCell2Wall(ARGS ... args)
        {
          cell2Wall = Node2NodeForce(std::forward<ARGS>(args)...);
        }

      protected:
        //! All lattice information and then some
        geometry::LatticeData &latticeData;
        //! Contains all cells
        CellContainer cells;
        //! Divide and conquer object
        DivideConquerCells cellDnC;
        //! Divide and conquer object
        DivideConquer<WallNode> wallDnC;
        //! A work array with forces/positions
        std::vector<LatticePosition> work;
        //! This function is called every lb turn
        //! It should insert cells using the call back passed to it.
        std::function<void(CellInserter const&)> cellInsertionCallBack;
        //! Observers to be notified when cell positions change
        std::vector<CellChangeListener> cellChangeListeners;
        //! Remove cells if they reach these outlets
        std::vector<FlowExtension> outlets;
        //! Interaction terms between cells
        Node2NodeForce cell2Cell;
        //! Interaction terms between cells
        Node2NodeForce cell2Wall;
    };

    template<class TRAITS>
    void CellArmy<TRAITS>::Fluid2CellInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Fluid -> cell interations");
      std::vector<LatticePosition> & positions = work;
      LatticePosition const origin(0, 0, 0);

      CellContainer::const_iterator i_first = cells.begin();
      CellContainer::const_iterator const i_end = cells.end();
      for (; i_first != i_end; ++i_first)
      {
        positions.resize( (*i_first)->GetVertices().size());
        std::fill(positions.begin(), positions.end(), origin);
        velocitiesOnMesh<Kernel, Stencil>(*i_first, latticeData, positions);
        (*i_first)->operator+=(positions);
      }
      // Positions have changed: update Divide and Conquer stuff
      cellDnC.update();
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::Cell2FluidInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Cell -> fluid interations");
      latticeData.ResetForces();

      for (auto const &cell: cells)
      {
        forcesOnGrid<typename Kernel::LatticeType, Stencil>(cell, work, latticeData);
      }
      addCell2CellInteractions<Stencil>(cellDnC, cell2Cell, latticeData);
      addCell2WallInteractions<Stencil>(cellDnC, wallDnC, cell2Wall, latticeData);
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::CellRemoval()
    {
      auto i_first = cells.cbegin();
      auto const i_end = cells.cend();
      while (i_first != i_end)
      {
        auto const barycenter = (*i_first)->GetBarycenter();
        auto checkCell = [&barycenter](FlowExtension const &flow)
        {
          return contains(flow, barycenter);
        };
        // save current iterator and increment before potential removal.
        // removing the cell from the set should invalidate only the relevant iterator.
        auto const i_current = i_first;
        ++i_first;
        if (std::find_if(outlets.begin(), outlets.end(), checkCell) != outlets.end())
        {
          log::Logger::Log<log::Debug, log::OnePerCore>(
              "Removing cell at (%f, %f, %f)",
              barycenter.x, barycenter.y, barycenter.z
          );
          cellDnC.remove(*i_current);
          cells.erase(i_current);
        }
      }
    }
  }
}

#endif
