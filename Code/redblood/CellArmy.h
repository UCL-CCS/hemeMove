// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLARMY_H
#define HEMELB_REDBLOOD_CELLARMY_H

#include <algorithm>
#include <vector>
#include <memory>
#include <iomanip>

#include <boost/uuid/uuid_io.hpp>

#include "Exception.h"
#include "geometry/Domain.h"
#include "redblood/Cell.h"
#include "redblood/CellCell.h"
#include "redblood/WallCellPairIterator.h"
#include "redblood/GridAndCell.h"
#include "redblood/FlowExtension.h"
#include "redblood/types.h"
#include "redblood/parallel/SpreadForces.h"
#include "redblood/parallel/IntegrateVelocities.h"
#include "reporting/Timers.h"

namespace hemelb
{
  namespace redblood
  {
    //! \brief Federates the cells together so we can apply ops simultaneously
    //! \tparam TRAITS holds type of kernel and stencil
    template<class TRAITS> class CellArmy
    {
      public:
        using Lattice = typename TRAITS::Lattice;
        using Kernel = typename TRAITS::Kernel;
        using Stencil = typename TRAITS::Stencil;

        CellArmy(geometry::FieldData &latDat, CellContainer const &cells,
                 std::shared_ptr<TemplateCellContainer> cellTemplates,
                 hemelb::reporting::Timers &timings, LatticeDistance boxsize = 10.0,
                 Node2NodeForce const &cell2Cell = { 0e0, 1e0, 2 },
                 Node2NodeForce const &cell2Wall = { 0e0, 1e0, 2 },
                 net::MpiCommunicator const &worldCommunicator = net::MpiCommunicator::World()) :
            fieldData(latDat), cells(cells), cellDnC(cells, boxsize, cell2Cell.cutoff + 1e-6),
                wallDnC(createWallNodeDnC<Lattice>(latDat.GetDomain(), boxsize, cell2Wall.cutoff + 1e-6)),
                cell2Cell(cell2Cell), cell2Wall(cell2Wall),
                cellTemplates(cellTemplates), timings(timings),
                neighbourDependenciesGraph(parallel::CreateGraphComm(worldCommunicator,
                                                                     latDat.GetDomain(),
                                                                     cellTemplates,
                                                                     timings)),
                exchangeCells(neighbourDependenciesGraph),
                velocityIntegrator(neighbourDependenciesGraph),
                forceSpreader(neighbourDependenciesGraph),
                globalCoordsToProcMap(parallel::ComputeGlobalCoordsToProcMap(neighbourDependenciesGraph, fieldData.GetDomain())),
                nodeDistributions(parallel::nodeDistributions(globalCoordsToProcMap, cells))
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
        parallel::LentCells const & GetLentCells() const
        {
          return lentCells;
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

        std::function<void(CellInserter const&)> GetCellInsertion() const
        {
          return cellInsertionCallBack;
        }

        //! Calls cell insertion
        void CallCellInsertion()
        {
          timings.cellInsertion().Start();
          if (cellInsertionCallBack)
          {
            auto callback = [this](CellContainer::value_type cell)
            {
              this->AddCell(cell);
            };
            cellInsertionCallBack(callback);
          }
          timings.cellInsertion().Stop();
        }

        //! Adds a cell change listener to be notified when cell positions change
        void AddCellChangeListener(CellChangeListener const & ccl)
        {
          cellChangeListeners.push_back(ccl);
        }

        //! Invokes the callback function to output cell positions
        void NotifyCellChangeListeners()
        {
          timings.cellListeners().Start();
          for (CellChangeListener ccl : cellChangeListeners)
          {
            ccl(cells);
          }
          timings.cellListeners().Stop();
        }

        //! Sets outlets within which cells disappear
        void SetOutlets(std::vector<FlowExtension> olets)
        {
          outlets = std::move(olets);
        }

        //! Remove cells if they have reached outlets
        void CellRemoval();

        //! Adds input cell to simulation
        void AddCell(CellContainer::value_type cell);

        //! \brief Sets cell to cell interaction forces
        //! \details Forwards arguments to Node2NodeForce constructor.
        template<class ... ARGS> void SetCell2Cell(ARGS && ... args)
        {
          cell2Cell = Node2NodeForce(std::forward<ARGS>(args)...);
          cellDnC.SetBoxSizeAndHalo(cellDnC.GetBoxSize(), cell2Cell.cutoff + 1e-6);
        }
        //! \brief Sets cell to cell interaction forces
        //! \details Forwards arguments to Node2NodeForce constructor.
        template<class ... ARGS> void SetCell2Wall(ARGS && ... args)
        {
          cell2Wall = Node2NodeForce(std::forward<ARGS>(args)...);
          wallDnC = createWallNodeDnC<Lattice>(wallDnC,
                                               wallDnC.GetBoxSize(),
                                               cell2Wall.cutoff + 1e-6);
        }

      protected:
        //! All lattice information and then some
        geometry::FieldData &fieldData;
        //! Contains all cells
        CellContainer cells;
        //! Cells lent to this process
        parallel::LentCells lentCells;
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
        //! Container with the templates used to create the RBC meshes
        std::shared_ptr<TemplateCellContainer> cellTemplates;
        //! Timers object used to time different code sections
        hemelb::reporting::Timers &timings;
        //! Communicator defining the data dependencies between processors for spreading/interpolation
        net::MpiCommunicator neighbourDependenciesGraph;
        //! Exchange cells object
        parallel::ExchangeCells exchangeCells;
        //! Velocity integrator object
        parallel::IntegrateVelocities velocityIntegrator;
        //! Force spreader object
        parallel::SpreadForces forceSpreader;
        //! Map allowing to look up which of the neighbours in neighbourDependenciesGraph owns a given lattice site (given by global coordinates)
        parallel::GlobalCoordsToProcMap globalCoordsToProcMap;
        //! Object describing how the cells affect different subdomains
        parallel::NodeDistributions nodeDistributions;

    };

    template<class TRAITS>
    void CellArmy<TRAITS>::Fluid2CellInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Fluid -> cell interations");

      timings.exchangeCells().Start();
      auto ownership = [this](CellContainer::value_type cell)
      {
        auto owner = nodeDistributions.at(cell->GetTag()).DominantAffectedProc();
        if (owner == -1)
        {
          throw Exception() << "Process " << fieldData.GetDomain().GetCommunicator().Rank()
                            << " cannot determine the owner of cell " << cell->GetTag()
                            << " with barycentre " << cell->GetBarycentre();
        }
        return owner;
      };
      exchangeCells.PostCellMessageLength(nodeDistributions, cells, ownership);
      exchangeCells.PostCells(nodeDistributions, cells, ownership);
      auto const distCells = exchangeCells.ReceiveCells(cellTemplates);
      exchangeCells.Update(cells, distCells);
      exchangeCells.Update(nodeDistributions,
                           distCells,
                           parallel::details::AssessMPIFunction<Stencil>(globalCoordsToProcMap));
      timings.exchangeCells().Stop();

      // Actually perform velocity integration
      timings.computeAndPostVelocities().Start();
      velocityIntegrator.PostMessageLength(std::get<2>(distCells));
      velocityIntegrator.ComputeLocalVelocitiesAndUpdatePositions<TRAITS>(fieldData, cells);
      velocityIntegrator.PostVelocities<TRAITS>(fieldData, std::get<2>(distCells));
      timings.computeAndPostVelocities().Stop();

      timings.receiveVelocitiesAndUpdate().Start();
      velocityIntegrator.UpdatePositionsNonLocal(nodeDistributions, cells);
      timings.receiveVelocitiesAndUpdate().Stop();

      // Positions have changed: update node distributions
      timings.computeNodeDistributions().Start();
      for (auto cell : cells)
      {
        try
        {
          nodeDistributions.at(cell->GetTag()).template Reindex<Stencil>(globalCoordsToProcMap, cell);
        }
        catch (std::exception const &e)
        {
          // If the simulation went unstable, some vertices may have been advected out of the domain.
          // This will be picked up in Reindex when figuring out vertex ownership. The code will throw
          // to give us a chance to e.g. write all the cells to disk for debugging purposes

          std::stringstream message;
          message << "Mesh vertex belongs to cell " << cell->GetTag();
          log::Logger::Log<log::Error, log::OnePerCore>(message.str());

          throw;
        }
      }
      timings.computeNodeDistributions().Stop();

      // Positions have changed: update Divide and Conquer stuff
      log::Logger::Log<log::Debug, log::OnePerCore>("Number of lent cells: %i",
                                                    std::get<2>(distCells).size());
      timings.updateDNC().Start();
      cellDnC.update(distCells);
      timings.updateDNC().Stop();
      lentCells = std::move(std::get<2>(distCells));
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::Cell2FluidInteractions()
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Cell -> fluid interations");
      fieldData.ResetForces();

      timings.computeAndPostForces().Start();
      forceSpreader.PostMessageLength(nodeDistributions, cells);
      forceSpreader.ComputeForces(cells);
      forceSpreader.PostForcesAndNodes(nodeDistributions, cells);
      forceSpreader.SpreadLocalForces<TRAITS>(fieldData, cells);
      timings.computeAndPostForces().Stop();

      timings.receiveForcesAndUpdate().Start();
      forceSpreader.SpreadNonLocalForces<TRAITS>(fieldData);
      timings.receiveForcesAndUpdate().Stop();

      //! @todo Any changes required for these lines when running in parallel?
      timings.updateCellAndWallInteractions().Start();
      addCell2CellInteractions<Stencil>(cellDnC, cell2Cell, fieldData);
      addCell2WallInteractions<Stencil>(cellDnC, wallDnC, cell2Wall, fieldData);
      timings.updateCellAndWallInteractions().Stop();
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::CellRemoval()
    {
      timings.cellRemoval().Start();
      auto i_first = cells.cbegin();
      auto const i_end = cells.cend();
      while (i_first != i_end)
      {
        auto const barycentre = (*i_first)->GetBarycentre();
        auto checkCell = [&barycentre](FlowExtension const &flow)
        {
          return contains(flow, barycentre);
        };
        // save current iterator and increment before potential removal.
        // removing the cell from the set should invalidate only the relevant iterator.
        auto const i_current = i_first;
        ++i_first;
        if (std::find_if(outlets.begin(), outlets.end(), checkCell) != outlets.end())
        {
          std::stringstream message;
          message << "Removing cell "<< (*i_current)->GetTag() << " at " << barycentre;
          log::Logger::Log<log::Info, log::OnePerCore>(message.str());

          cellDnC.remove(*i_current);
          auto const numErased = nodeDistributions.erase((*i_current)->GetTag());
          cells.erase(i_current);
          assert(numErased == 1);
        }
      }
      timings.cellRemoval().Stop();
    }

    template<class TRAITS>
    void CellArmy<TRAITS>::AddCell(CellContainer::value_type cell)
    {
      auto const barycentre = cell->GetBarycentre();

      //! @todo: #623 AddCell should only be called if the subdomain contains the relevant RBC inlet
      // TODO: #759 truncation of barycentre
      auto const iter = globalCoordsToProcMap.find(Vec16{barycentre});
      bool insertAtThisRank = (iter != globalCoordsToProcMap.end()) && (iter->second == neighbourDependenciesGraph.Rank());
      if (insertAtThisRank)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Adding cell at (%f, %f, %f)",
            barycentre.x(),
            barycentre.y(),
            barycentre.z());
        cellDnC.insert(cell);
        cells.insert(cell);

        nodeDistributions.emplace(std::piecewise_construct,
            std::forward_as_tuple(cell->GetTag()),
            std::forward_as_tuple(parallel::details::AssessMPIFunction<
              Stencil>(globalCoordsToProcMap),
              cell));
        log::Logger::Log<log::Info, log::OnePerCore>("Cell has %i edge nodes",
          nodeDistributions.find(cell->GetTag())->second.BoundaryIndices().size());
      }

#ifndef NDEBUG
      // Check that one and only one process inserted the cell
      unsigned numCellsAdded = neighbourDependenciesGraph.AllReduce((unsigned) insertAtThisRank, MPI_SUM);
      if (numCellsAdded != 1)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Failed to add cell at (%f, %f, %f). It was added %d times.",
            barycentre.x(),
            barycentre.y(),
            barycentre.z(),
            numCellsAdded);

        hemelb::net::MpiEnvironment::Abort(-1);
      }
#endif

    }
  }
}

#endif
