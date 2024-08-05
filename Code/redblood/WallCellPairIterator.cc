// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <vector>
#include "redblood/WallCellPairIterator.h"

namespace hemelb::redblood
{
    bool WallCellPairIterator::operator++()
    {
      do
      {
        if (firstWallNode == wallNodes.end())
        {
          return false;
        }
        if (firstCellNode != lastCellNode)
        {
          ++firstCellNode;
        }
        else if (++box_iterator)
        {
          std::tie(firstCellNode, lastCellNode) = cellNodes(firstWallNode->second.node
              + *box_iterator);
        }
        else if (++firstWallNode != wallNodes.end())
        {
          box_iterator = BorderBoxIterator(firstWallNode->second.nearBorder);
          // should always at least include cell of firstWallNode
          assert(static_cast<bool>(box_iterator));
          std::tie(firstCellNode, lastCellNode) = cellNodes(firstWallNode->second.node
              + *box_iterator);
        }
        else
        {
          // reached end of loop
          return false;
        }
      }
      while (not isValid());
      return true;
    }

    bool WallCellPairIterator::isValid() const
    {
      return firstCellNode != lastCellNode and box_iterator
          and (*firstCellNode - firstWallNode->second.node).GetMagnitude() < cutoff;
    }

    WallCellPairIterationRange iterate(
            DivideConquerCells const& cellDnC,
            DivideConquer<WallNode> const& wallDnC,
            LatticeDistance const & cutoff)
    {
        return {&cellDnC, &wallDnC, cutoff};
    }

    WallCellPairIterator WallCellPairIterationRange::begin() {
        return {*cellDnC, *wallDnC, cutoff, WallCellPairIterator::Begin()};
    }
    WallCellPairIterator WallCellPairIterationRange::end() {
        return {*cellDnC, *wallDnC, cutoff, WallCellPairIterator::End()};
    }

    WallCellPairIterator::WallCellPairIterator(
        DivideConquerCells const& cellNodes, DivideConquer<WallNode> const &wallNodes,
        LatticeDistance cutoff, DivideConquerCells::const_iterator const &firstCellNode,
        DivideConquerCells::const_iterator const &lastCellNode,
        DivideConquer<WallNode>::const_iterator const &firstWallNode) :
        cellNodes(cellNodes), wallNodes(wallNodes), firstCellNode(firstCellNode),
            lastCellNode(lastCellNode), firstWallNode(firstWallNode), box_iterator(0),
            cutoff(cutoff)
    {
      if (firstWallNode == wallNodes.end())
      {
        return;
      }
      try
      {
        box_iterator = BorderBoxIterator(firstWallNode->second.nearBorder);
        assert(static_cast<bool>(box_iterator));
        if (firstCellNode == lastCellNode)
        {
          std::tie(this->firstCellNode, this->lastCellNode) = cellNodes(firstWallNode->second.node
              + *box_iterator);
        }
        if (not isValid())
        {
          operator++();
        }
      }
      catch (...)
      {
        this->firstWallNode = wallNodes.end();
        throw;
      }
    }
}
