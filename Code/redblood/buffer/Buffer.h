// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_BUFFER_BUFFER_H
#define HEMELB_REDBLOOD_BUFFER_BUFFER_H

#include <utility>

#include "redblood/Cell.h"
#include "redblood/FlowExtension.h"
#include "redblood/types.h"
#include "units.h"

namespace hemelb::redblood::buffer
{
      //! Container of virtual cells
      class Buffer
      {
        public:
          //! Signature of a function returning a new cell at each call
          using CellDistributionFunction = std::function<CellContainer::value_type ()>;
          //! Constructs buffer from geometry and cells
          //! \param[in] cyl: Geometric description of the buffer
          //!     The origin of the cylinder is the drop point in the coordinates of LB. The normal
          //!     should point from the region nearest to drop point to furthest from drop point. It
          //!     is co-linear with the axis of the cylinder.
          //! \param[in] cells: Virtuals cells at input.
          Buffer(std::shared_ptr<Cylinder> cyl, CellContainer cells = CellContainer()) :
              geometry(std::move(cyl)), virtuals(std::move(cells)),
              offset(0), interactionRadius(0), numberOfRequests(0)
          {
          }
          //! Constructs buffer from geometry only
          Buffer(Cylinder const & cyl, CellContainer cells = CellContainer()) :
              geometry(std::make_shared<Cylinder>(cyl)), virtuals(std::move(cells)),
              offset(0), interactionRadius(0), numberOfRequests(0)
          {
          }
          //! Destroys buffer
          virtual ~Buffer() = default;

          /**
           * Cell insertion callback called on each step of the simulation.  Cells
           * are inserted into the inlets while condition() evaluates to true.
           * Multiple cells are potentially inserted on each iteration.
           *
           * @param insertFn the function to insert a new cell into the simulation
           *
           * @see hemelb::redblood::CellArmy::SetCellInsertion
           * @see hemelb::redblood::CellArmy::CallCellInsertion
           */
          void operator()(CellInserter const& insertFn);

          //! Request buffer to drop n new cells
          void requestNewCells(site_t n = 1)
          {
            numberOfRequests += n;
          }
          //! Number of cells requested for drop off
          site_t NumberOfRequests() const
          {
            return numberOfRequests;
          }
          void SetCellDistributionFunction(CellDistributionFunction const &func)
          {
            getNewVirtualCell = func;
          }

        protected:
          //! Minimum number of cells to have at the ready
          const site_t minCells = 3;
          //! Geometry of the buffer
          std::shared_ptr<Cylinder> geometry;
          //! Container of virtual cells
          CellContainer virtuals;
          //! Offset between LB and buffer coordinate systems.
          LatticeDistance offset;
          //! Last dropped
          CellContainer::value_type justDropped;
          //! Last dropped cell's position
          LatticeDistance lastZ;
          //! Cell interaction distance
          LatticeDistance interactionRadius;
          //! Number of requests cells for drop-off
          site_t numberOfRequests;
          //! when called, returns a new cell to add to the buffer
          CellDistributionFunction getNewVirtualCell;

          //! Inserts a cell. Does not check whether it is inside the geometry of the buffer, nor
          //! whether it overlaps with other cells.
          auto insert(const CellContainer::value_type& cell)
          {
            return virtuals.insert(cell);
          }
          //! Updates offset between LB and offset
          void updateOffset();
          //! Drops the next nearest cell
          CellContainer::value_type drop();
          //! True if next cell can be dropped
          bool isDroppablePosition(CellContainer::value_type const &candidate) const
          {
            return isDroppablePosition(candidate->GetBarycentre());
          }
          //! True if position corresponds to that of a droppable cell
          //! This function works in the buffer's cartesian coordinates. It will take care of adding
          //! the offset.
          bool isDroppablePosition(LatticePosition const &) const;
          //! Fills buffer with new cells
          //! A minimum of n cells are added. Then further cells are added until the cylinder is
          //! filled. Includes a interaction radius buffer, so that if a cell is dropped, it will
          //! have virtual neighbors that will keep it from backtracking.
          void fillBuffer(site_t n);
          //! Returns next cell to drop
          CellContainer::value_type nearestCell() const;
          //! Returns next cell to drop
          CellContainer::value_type furthestCell() const;
      };

} // namespace hemelb::redblood
#endif
