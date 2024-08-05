// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_WALLCELLPAIRITERATOR_H
#define HEMELB_REDBLOOD_WALLCELLPAIRITERATOR_H

#include <cassert>
#include <memory>
#include <tuple>
#include <iterator>
#include "units.h"
#include "redblood/DivideConquer.h"
#include "redblood/CellCell.h"
#include "redblood/Borders.h"
#include "geometry/Domain.h"

namespace hemelb::redblood
{
    //! References a node of a mesh in the divide-and-conquer box
    class WallNode
    {
    public:
        //! Index of node in mesh
        LatticePosition node;
        //! Whether the node is near the border of the cube;
        size_t nearBorder;
    };

    //! \brief Creates divide and conquer box of wall nodes
    //! \brief The items of the divide and conquer box contain the location of the wall-node (not
    //! the location of the associated fluid site) as well as how close the node is to neighboring
    //! boxes.
    template<class LATTICE>
    DivideConquer<WallNode> createWallNodeDnC(geometry::Domain const&domain,
                                              LatticeDistance boxSize,
                                              LatticeDistance interactionDistance);
    //! \brief Creates divide and conquer box of wall nodes
    //! \brief This version is really to update the size or halo of the DnC
    template<class LATTICE>
    DivideConquer<WallNode> createWallNodeDnC(std::vector<LatticePosition> const&,
                                              LatticeDistance boxSize,
                                              LatticeDistance interactionDistance);
    //! \brief Creates divide and conquer box of wall nodes
    //! \brief This version is really to update the size or halo of the DnC
    template<class LATTICE>
    DivideConquer<WallNode> createWallNodeDnC(DivideConquer<WallNode> const &,
                                              LatticeDistance boxSize,
                                              LatticeDistance interactionDistance);

    //! Iterates over pairs of wall-node, cell-node which are within a given distance
    class WallCellPairIterator
    {
    public:
        //! Underlying type of the dereferenced object
        struct value_type
        {
            LatticePosition const &cellNode;
            LatticePosition const &wallNode;
        };
        using pointer = std::unique_ptr<WallCellPairIterator::value_type>;
        using reference = WallCellPairIterator::value_type const&;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;
        //! Tag to initialize a begining-of-range iterator
        struct Begin
        {
        };
        //! Tag to initialize an end-of-range iterator
        struct End
        {
        };
        //! Creates an iterator of Wall and Cell nodes
        WallCellPairIterator(DivideConquerCells const& cellNodes,
                             DivideConquer<WallNode> const &wallNodes, LatticeDistance cutoff,
                             WallCellPairIterator::Begin const &) :
                WallCellPairIterator(cellNodes,
                                     wallNodes,
                                     cutoff,
                                     cellNodes.end(),
                                     cellNodes.end(),
                                     wallNodes.begin())
        {
        }
        WallCellPairIterator(DivideConquerCells const& cellNodes,
                             DivideConquer<WallNode> const &wallNodes, LatticeDistance cutoff) :
                WallCellPairIterator(cellNodes, wallNodes, cutoff, Begin())
        {
        }
        //! Creates an iterator of Wall and Cell nodes
        WallCellPairIterator(DivideConquerCells const& cellNodes,
                             DivideConquer<WallNode> const &wallNodes, LatticeDistance cutoff,
                             WallCellPairIterator::End const &) :
                WallCellPairIterator(cellNodes,
                                     wallNodes,
                                     cutoff,
                                     cellNodes.end(),
                                     cellNodes.end(),
                                     wallNodes.end())
        {
        }
        //! Positions of the wall and cell nodes
        value_type operator*() const
        {
            HASSERT(static_cast<bool>(*this));
            return {*firstCellNode, firstWallNode->second.node};
        }
        //! Do not use. Only exists to model an InputIterator
        pointer operator->() const
        {
            HASSERT(static_cast<bool>(*this));
            value_type tmp{*firstCellNode, firstWallNode->second.node};
            return std::make_unique<value_type>(tmp);
        }

        bool operator++();
        WallCellPairIterator operator++(int)
        {
            HASSERT(static_cast<bool>(*this));
            WallCellPairIterator result(cellNodes,
                                        wallNodes,
                                        cutoff,
                                        firstCellNode,
                                        lastCellNode,
                                        firstWallNode);
            operator++();
            return result;
        }

        operator bool() const
        {
            return firstWallNode != wallNodes.end() and isValid();
        }

    protected:
        //! Reference to node vertices via divide and conquer boxes
        DivideConquerCells const& cellNodes;
        //! wall node iterator
        DivideConquer<WallNode> const& wallNodes;
        //! cell node iterator
        DivideConquerCells::const_iterator firstCellNode;
        //! cell node iterator
        DivideConquerCells::const_iterator lastCellNode;
        //! wall node iterator
        DivideConquer<WallNode>::const_iterator firstWallNode;
        //! iterates over difference boxes around current wall node
        BorderBoxIterator box_iterator;
        //! Maximum interaction distance
        LatticeDistance cutoff;

        //! Whether current pair is within range
        bool isValid() const;

        //! Creates an iterator of Wall and Cell nodes
        WallCellPairIterator(DivideConquerCells const& cellNodes,
                             DivideConquer<WallNode> const &wallNodes, LatticeDistance cutoff,
                             DivideConquerCells::const_iterator const &firstCellNode,
                             DivideConquerCells::const_iterator const &lastCellNode,
                             DivideConquer<WallNode>::const_iterator const &firstWallNode);
    };

    template<class LATTICE>
    DivideConquer<WallNode> createWallNodeDnC(geometry::Domain const&domain,
                                              LatticeDistance boxSize,
                                              LatticeDistance interactionDistance)
    {
        DivideConquer<WallNode> result(boxSize);
        for (site_t i = 0; i < domain.GetLocalFluidSiteCount(); ++i) {
            auto const site = domain.GetSite(i);
            if (!site.IsWall())
                continue;
            for (Direction ii = 1; ii < LATTICE::NUMVECTORS; ++ii) {
                if (!site.HasWall(ii))
                    continue;

                LatticeDistance const distance = site.GetWallDistance<LATTICE>(ii);

                // Direction of streaming from wall to this site
                auto const& direction = LATTICE::CD[ii];
                auto const wallnode = LatticePosition(site.GetGlobalSiteCoords()) + direction.GetNormalised() * distance;
                auto const nearness = figureNearness(result, wallnode, interactionDistance);
                result.insert(wallnode, { wallnode, nearness });
            }
        }
        return result;
    }

    template<class LATTICE>
    DivideConquer<WallNode> createWallNodeDnC(std::vector<LatticePosition> const& nodes,
                                              LatticeDistance boxSize,
                                              LatticeDistance interactionDistance)
    {
        DivideConquer<WallNode> result(boxSize);
        for (auto&& node: nodes) {
            auto const nearness = figureNearness(result, node, interactionDistance);
            result.insert(node, {node, nearness});
        }
        return result;
    }

    template<class LATTICE>
    DivideConquer<WallNode> createWallNodeDnC(DivideConquer<WallNode> const& nodes,
                                              LatticeDistance boxSize,
                                              LatticeDistance interactionDistance)
    {
        DivideConquer<WallNode> result(boxSize);
        for (auto&& node: nodes) {
            auto const nearness = figureNearness(result, node.second.node, interactionDistance);
            result.insert(node.second.node, { node.second.node, nearness });
        }
        return result;
    }

    //! \brief Creates an object that std::begin and std::end can act on
    //! \details This makes it possible to use ranged for loops as well as iterate in more
    //! tradditional manner:
    //! \code{.cpp}
    //!   auto const last = std::end(iterate(cellDnC, wallDnC, interactionDistance));
    //!   auto first = std::begin(iterate(cellDnC, wallDnC, interactionDistance));
    //!   for(; first != last; ++first)
    //!     std::cout << first->wallNode << " " << first->cellNode << "\n";
    //| \endcode
    //! or
    //! \code{.cpp}
    //!   for(auto const nodes: iterate(cellDnC, wallDnC, interactionDistance))
    //!     std::cout << nodes.wallNode << " " << nodes.cellNode << "\n";
    //| \endcode
    struct WallCellPairIterationRange {
        DivideConquerCells const* cellDnC;
        DivideConquer<WallNode> const* wallDnC;
        LatticeDistance cutoff;

        WallCellPairIterator begin();
        WallCellPairIterator end();
    };

    WallCellPairIterationRange iterate(
            DivideConquerCells const& cellDnC,
            DivideConquer<WallNode> const& wallDnC,
            LatticeDistance const & cutoff);

    //! \brief Computes cell <-> cell interactions and spread to grid
    //! \details Given a partition of the cells' nodes and node <-> node interaction
    //! functional, computes the short-range that can occur between cells that are
    //! too close to one another. The interaction forces are computed and spread to
    //! the lattice.
    template<class STENCIL>
    void addCell2WallInteractions(DivideConquerCells const &cellDnC,
                                  DivideConquer<WallNode> const &wallDnC,
                                  Node2NodeForce const &functional,
                                  geometry::FieldData &latticeData);

    // implementation must happen after begin and end have been declared
    template<class STENCIL>
    void addCell2WallInteractions(DivideConquerCells const &cellDnC,
                                  DivideConquer<WallNode> const &wallDnC,
                                  Node2NodeForce const &functional,
                                  geometry::FieldData &latticeData)
    {
        for (WallCellPairIterator iter { cellDnC,
                                         wallDnC,
                                         functional.cutoff,
                                         WallCellPairIterator::Begin() }; iter; ++iter)
        {
            auto const force = functional(iter->cellNode, iter->wallNode);
            spreadForce<STENCIL>(iter->cellNode, latticeData, force);
        }
    }
}

#endif
