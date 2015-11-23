//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//
#ifndef HEMELB_REDBLOODE_PARALLELIZATION_CELL_PARALLELIZATION
#define HEMELB_REDBLOODE_PARALLELIZATION_CELL_PARALLELIZATION

#include <boost/uuid/uuid.hpp>
#include <map>

#include "redblood/parallel/NodeCharacterizer.h"
#include "redblood/Cell.h"
#include "net/MpiCommunicator.h"
#include "net/INeighborAllToAll.h"
#include "net/INeighborAllToAllV.h"


namespace hemelb
{
  namespace redblood
  {
    namespace parallel
    {
      class CellParallelization
      {
        public:
          //! Type of the object holding distributions
          typedef std::map<boost::uuids::uuid, NodeCharacterizer> NodeDistributions;
          //! Container holding cells lent by other processes
          typedef std::map<proc_t, CellContainer> LentCells;
          //! Creates a cell-parallelization object
          CellParallelization(net::MpiCommunicator const &comm)
            : comm(comm.Duplicate())
          {
          }
          //! Creates a cell-parallelization object
          CellParallelization(net::MpiCommunicator &&comm)
            : comm(std::move(comm))
          {
          }
          //! Adds an owned cell
          void AddCell(CellContainer::const_reference cell)
          {
            owned.insert(cell);
          }

        protected:
          //! Graph communicator defining neighberhood over which cells can be owned
          net::MpiCommunicator comm;
          //! Cells owned by this process
          CellContainer owned;
          //! Distribution of the cells owned by this process
          NodeDistributions distributions;
          //! Cells lent by other proces
          LentCells lent;
      };

      class ExchangeCells
      {
        public:
          //! Type of the object holding distributions
          typedef CellParallelization::NodeDistributions NodeDistributions;
          //! Container holding cells lent by other processes
          typedef CellParallelization::LentCells LentCells;
          //! Function that can figure out who should own a cell
          typedef std::function<int(CellContainer::const_reference)> Ownership;
          //! Result of the whole messaging mess
          typedef std::tuple<CellContainer, CellContainer, LentCells> ChangedCells;

          ExchangeCells(net::MpiCommunicator const &graphComm)
            : cellCount(graphComm), totalNodeCount(graphComm), nameLengths(graphComm),
              templateNames(graphComm), ownerIDs(graphComm), nodeCount(graphComm),
              cellUUIDs(graphComm), cellScales(graphComm), nodePositions(graphComm)
          {
          }
          //! \brief Computes and posts length of message when sending cells
          //! \param[in] distributions: Node distributions of the cells owned by this process
          //! \param[in] owned: Cells currently owned by this process
          //! \param[in] ownership a function to ascertain ownership. It should return the rank of
          //! the owning process in the *world* communicator, according to the position in the cell.
          virtual void PostCellMessageLength(
              NodeDistributions const& distributions, CellContainer const &owned,
              Ownership const & ownership);
          //! \brief Computes and posts length of message when sending cells
          //! \param[in] distributions: Node distributions of the cells owned by this process
          //! \param[in] owned: Cells currently owned by this process
          //! \param[in] ownership id of the process owning the cell, corresponding to the rank in
          //! the graph communicator used to build this object.
          virtual void PostCellMessageLength(
              NodeDistributions const& distributions, CellContainer const &owned,
              std::map<boost::uuids::uuid, proc_t> const & ownership);
          //! \brief Post all owned cells and preps for receiving lent cells
          //! \param[in] distributions tells us for each proc the list of nodes it requires
          //! \param[in] cells a container of cells owned and managed by this process
          //! \param[in] ownership a function to ascertain ownership. It should return the rank of
          //! the owning process in the *world* communicator, according to the position in the cell.
          virtual void PostCells(
              NodeDistributions const &distributions, CellContainer const &cells,
              Ownership const & ownership);
          //! \brief Post all owned cells and preps for receiving lent cells
          //! \param[in] distributions tells us for each proc the list of nodes it requires
          //! \param[in] cells a container of cells owned and managed by this process
          //! \param[in] ownership id of the process owning the cell, corresponding to the rank in
          //! the graph communicator used to build this object.
          virtual void PostCells(
              NodeDistributions const &distributions, CellContainer const &cells,
              std::map<boost::uuids::uuid, proc_t> const & ownership);
          //! Receives messages, reconstructs cells
          //! \return a 3-tuple with the newly owned cells, the disowned cells, and the lent cells
          ChangedCells ReceiveCells(
              std::shared_ptr<TemplateCellContainer const> const &templateCells)
          {
            return ReceiveCells(*templateCells);
          }
          virtual ChangedCells ReceiveCells(TemplateCellContainer const &templateCells);

          //! Adds new cells and removes old ones
          static void Update(CellContainer &owned, ChangedCells const & changes);

          //! Updates node-distributions to have only owned cells
          static void Update(
              NodeDistributions &distributions, ChangedCells const & changes,
              NodeCharacterizer::AssessNodeRange const &assessor);
        protected:
          //! \brief Sends number of cells
          //! \details Using int because it meshes better with sending the number of nodes per cell.
          net::INeighborAllToAll<int> cellCount;
          //! Sends total number of nodes
          net::INeighborAllToAll<int> totalNodeCount;
          //! Send size of message containing template mesh names
          net::INeighborAllToAll<int> nameLengths;
          //! Names of the template meshes
          net::INeighborAllToAllV<char> templateNames;
          //! ID of owner process
          net::INeighborAllToAllV<int> ownerIDs;
          //! Sends number of nodes per cell
          net::INeighborAllToAllV<size_t> nodeCount;
          //! Cell uuids
          net::INeighborAllToAllV<unsigned char> cellUUIDs;
          //! Sends scales
          net::INeighborAllToAllV<LatticeDistance> cellScales;
          //! Sends positions
          net::INeighborAllToAllV<LatticeDistance> nodePositions;
          //! \brief Cell that are no longuer owned by this process
          //! \details Unlike formelyOwned, this keeps track of the whole cell
          CellContainer disowned;
          //! \brief Formely owned cells are lent back to this process by the neibhbor
          //! \details These cells are the same as the disowned cells. However, only part of the
          //! nodes kept: those that affect this process.
          LentCells formelyOwned;

          //! Number of nodes to send to each neighboring process
          void SetupLocalSendBuffers(
              NodeDistributions const &distributions, CellContainer const &cells,
              std::map<boost::uuids::uuid, proc_t> const & ownership);
          //! Adds data to local buffers for a cell that retains the same ownership
          void AddOwnedToLocalSendBuffers(
              int neighbor, int nth, int nVertices, proc_t owner,
              CellContainer::const_reference cell,
              NodeCharacterizer::Process2NodesMap::mapped_type const& indices);
          //! Adds data to local buffer for a cell that changes ownership
          void AddDisownedToLocalSendBuffers(
              int neighbor, int nth, CellContainer::const_reference cell);
          //! Adds data to local buffer for a cell that changes ownership
          void AddDisownedToLocalSendBuffers(
              int neighbor, int nth, CellContainer::const_reference cell,
              NodeCharacterizer::Process2NodesMap::mapped_type const& indices);
          //! Adds local data exept nodes
          void AddToLocalSendBuffersAllButNodes(
              int neighbor, int nth, proc_t ownerID, CellContainer::const_reference cell);
          //! Recreates a given cell from messages
          CellContainer::value_type RecreateLentCell(size_t i);
          //! Recreates a given cell from messages
          CellContainer::value_type RecreateOwnedCell(
              size_t i, TemplateCellContainer const &templateCells);
      };

    } /* parallel */
  } /* redblood */
} /* hemelb */
#endif
