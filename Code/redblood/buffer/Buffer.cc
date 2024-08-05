// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iterator>
#include "redblood/buffer/Buffer.h"
#include "redblood/types.h"
#include "Exception.h"

namespace hemelb::redblood::buffer
{
      namespace
      {
        LatticeDistance maxCellRadius(CellBase const& cell)
        {
          auto const barycentre = cell.GetBarycentre();
          auto const &vertices = cell.GetVertices();
          auto const first = vertices.begin();
          auto dist = [&barycentre](LatticePosition const &a, LatticePosition const &b)
          {
            return (a-barycentre).GetMagnitudeSquared() < (b-barycentre).GetMagnitudeSquared();
          };
          LatticePosition const &max = *std::max_element(first, vertices.end(), dist);
          return 4e0 * (max - barycentre).GetMagnitudeSquared();
        }

        template<class T_FUNC>
        CellContainer::value_type orderedCell(T_FUNC measure, CellContainer const &cells)
        {
          if (cells.empty())
            throw Exception() << "No cells in buffer";

          std::vector<LatticeDistance> dists(cells.size(), 0);
          std::transform(cells.begin(), cells.end(), dists.begin(), measure);
          auto const n = std::min_element(dists.begin(), dists.end()) - dists.begin();
          auto min_iter = cells.begin();
          std::advance(min_iter, n);
          return *min_iter;
        }
      }

      CellContainer::value_type Buffer::nearestCell() const
      {
        // normal should point from nearest to drop point to furthest from drop point
        // e.g. from origin to inside of the buffer
        // e.g. opposite to vascular system
        // The drop point is *always* (0, 0, 0) in the coordinate system of the virtual buffer.
        auto const &normal = geometry->normal;
        auto getdist = [&normal](CellContainer::value_type const& c)
        {
          return Dot(c->GetBarycentre(), normal);
        };
        return orderedCell(getdist, virtuals);
      }

      CellContainer::value_type Buffer::furthestCell() const
      {
        auto const &normal = geometry->normal;
        auto getdist = [&normal](const CellContainer::value_type& c)
        {
          return -Dot(c->GetBarycentre(), normal);
        };
        return orderedCell(getdist, virtuals);
      }

      CellContainer::value_type Buffer::drop()
      {
        justDropped = nearestCell();
        virtuals.erase(justDropped);

        lastZ = Dot(justDropped->GetBarycentre(), geometry->normal);
        *justDropped += geometry->origin + geometry->normal * offset;
        return justDropped;
      }

      bool Buffer::isDroppablePosition(LatticePosition const &position) const
      {
        // position with respect to cylinder updated by offset
        auto const z = Dot(position + geometry->normal * offset, geometry->normal);
        return z <= geometry->length;
      }

      void Buffer::fillBuffer(site_t n)
      {
        if (not getNewVirtualCell)
        {
          throw Exception() << "Function to fill virtual buffer with cells is not set";
        }
        // Function that inserts and returns new cell
        auto insertCell = [this]()
        {
          auto const a = this->getNewVirtualCell();
          return *(this->virtuals.insert(a).first);
        };
        // Make sure buffer is not empty
        if (virtuals.empty() and n <= 0)
        {
          n = 1;
        }
        // Add the requested cells first
        CellContainer::value_type lastCell;
        for (site_t i(0); i < n; ++i)
        {
          lastCell = insertCell();
        }
        // if no cells were added, then find the cell furthest from being dropped
        if (not lastCell)
        {
          lastCell = furthestCell();
        }
        // add cell until outside geometry, including interaction radius buffer.
        while (isDroppablePosition(lastCell->GetBarycentre() - geometry->normal * interactionRadius))
        {
          lastCell = insertCell();
        }
      }

      void Buffer::operator()(CellInserter const& insertFn)
      {
        // Add virtual cells, if necessary
        if (virtuals.size() < static_cast<std::size_t>(NumberOfRequests()))
        {
          fillBuffer(NumberOfRequests() - ssize(virtuals));
        }
        // Drop as many cells as possible
        for (; NumberOfRequests() > 0; --numberOfRequests)
        {
          // First update offsets
          updateOffset();
          if (not isDroppablePosition(nearestCell()))
          {
            break;
          }
          // Add to wherever
          insertFn(drop());
        }
        // Make sure buffer is still filled to prevent newly dropped cells from moving into space
        // that should be occupied by virtual cells.
        updateOffset();
        fillBuffer(0);
      }

      void Buffer::updateOffset()
      {
        // if no cell, then reset offset to zero
        if (virtuals.empty())
        {
          offset = 0;
          return;
        }

        // Makes sure there is an interaction distance, otherwise compute it
        if (interactionRadius <= 0e0)
        {
          interactionRadius = 1.25 * maxCellRadius(**virtuals.begin());
        }

        auto const normal = geometry->normal;
        auto const zCell = Dot(normal, nearestCell()->GetBarycentre());
        if (not justDropped)
        {
          offset = -zCell;
          return;
        }
        auto const zDropped = Dot(normal, justDropped->GetBarycentre() - geometry->origin) - offset;
        // dropped cell moved forward but distance with next cell is small
        // Then movement must the smallest of:
        // - distance moved by dropped cell over last LB iteration
        // - subject to the condition that the next cell does not move past the drop point.
        // The former will be negative. We've already checked that it is.
        if (zCell - interactionRadius < zDropped)
        {
          auto const delta = lastZ - zDropped;
          lastZ = zDropped;
          // Only move if dropped cell moved forward
          if (delta > 0)
          {
            offset -= std::min(delta, zCell);
          }
          return;
        }

        // Now we can move by as much as zCell, as long as we don't reach the interaction radius
        offset -= std::min(zCell, zCell - zDropped - interactionRadius);
      }
} // hemelb::redblood
