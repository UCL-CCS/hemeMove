//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARALLEL_CELLPARALLELIZATION_H
#define HEMELB_UNITTESTS_REDBLOOD_PARALLEL_CELLPARALLELIZATION_H

#include <cppunit/TestFixture.h>
#include <functional>

#include "unittests/redblood/Fixtures.h"
#include "net/MpiCommunicator.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "redblood/parallel/CellParallelization.h"
#include "util/Iterator.h"
#include <algorithm>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood_parallel
    {
      class ExchangeCells : public redblood::parallel::ExchangeCells
      {
        public:
          ExchangeCells(net::MpiCommunicator const &graphComm)
            : redblood::parallel::ExchangeCells(graphComm)
          {
          };

          std::vector<Length> const& GetSendLengths() const
          {
            return sendLengths;
          }
          std::vector<Length> const& GetReceiveLengths() const
          {
            return receiveLengths;
          }
          void WaitOutLengthRequest()
          {
            MPI_Status status;
            HEMELB_MPI_CALL(MPI_Wait, (&lengthRequest, &status));
          }
      };

      using namespace hemelb::redblood;
      using namespace hemelb::redblood::parallel;
      class CellParallelizationTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (CellParallelizationTests);
          CPPUNIT_TEST(testCellSwapGetLength);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();
          //! Test message sending number of cells and number of nodes
          void testCellSwapGetLength();

          //! Set of nodes affected by given proc
          std::set<proc_t> nodeLocation(LatticePosition const &node);

          //! Center of the region affected by a node
          LatticePosition GetCenter(proc_t rank) const;

          //! Creates map of distributed nodes
          CellParallelization::NodeDistributions GetNodeDistribution(
              CellContainer const &cells) const;

          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(
              LatticePosition const &pos, Dimensionless scale=1e0, unsigned int depth=0) const;
          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(Dimensionless scale=1e0, unsigned int depth=0) const
          {
            return GetCell(GetCenter(graph ? graph.Rank(): 0), scale, depth);
          }

        protected:
          LatticeDistance const radius = 5;
          net::MpiCommunicator graph;
      };

      void CellParallelizationTests::setUp()
      {
        // setups graph communicator
        auto world = net::MpiCommunicator::World();
        if(world.Size() >= 4)
        {
          std::vector<std::vector<int>> vertices{{1}, {0, 2, 3}, {1, 3}, {1, 2}};
          for(int i(4); i < world.Size(); ++i)
          {
            vertices.push_back(std::vector<int>{});
          }
          graph = std::move(world.Graph(vertices));
        }
      }

      //! Get cell centered at given position
      std::shared_ptr<Cell> CellParallelizationTests::GetCell(
          LatticePosition const &pos, Dimensionless scale, unsigned int depth) const
      {
        auto cell = std::make_shared<Cell>(icoSphere(depth));
        *cell += pos - cell->GetBarycenter();
        *cell *= scale;
        return cell;
      }

      CellParallelization::NodeDistributions CellParallelizationTests::GetNodeDistribution(
          CellContainer const & ownedCells) const
      {
        CellParallelization::NodeDistributions result;
        auto const assessor = std::bind(
            &CellParallelizationTests::nodeLocation, *this, std::placeholders::_1);
        for(auto const&cell: ownedCells)
        {
          result.emplace(
              std::piecewise_construct,
              std::forward_as_tuple(cell->GetTag()),
              std::forward_as_tuple(assessor, cell)
          );
        }
        return result;
      }

      LatticePosition CellParallelizationTests::GetCenter(proc_t rank) const
      {
        switch(rank)
        {
          case 0:
          case 1:
            return {0, 0, static_cast<double>(rank * 2) * radius};
          case 2:
            return {0, 2 * radius, (static_cast<double>(rank * 2 + 1)) * radius};
          case 3:
            return {0, 0, (static_cast<double>(rank * 2 - 2)) * radius};
          default:
            return {0, 0, static_cast<double>(rank) * radius * 10};
        }
      }

      std::set<proc_t> CellParallelizationTests::nodeLocation(LatticePosition const &node)
      {
        std::set<proc_t> result;
        for(proc_t i(0); i < graph.Size(); ++i)
        {
          if((node - GetCenter(i)).GetMagnitude() < radius)
          {
            result.insert(i);
          }
        }
        return result;
      }

      //! Check that cell send each other whole cells
      void CellParallelizationTests::testCellSwapGetLength()
      {
        if(not graph)
        {
          return;
        }
        size_t const sendto =
          graph.Rank() == 0 ? std::numeric_limits<size_t>::max():
          graph.Rank() == 1 ? 0:
          graph.Rank() == 2 ? 1:
          graph.Rank() == 3 ? 2: std::numeric_limits<size_t>::max();
        auto const center = GetCenter(sendto);
        CellContainer cells{ GetCell(center, 1e0, graph.Rank()) };
        if(graph.Rank() == 0)
        {
          cells.clear();
        }
        auto const dist = GetNodeDistribution(cells);

        ExchangeCells xc(graph);
        xc.PostCellMessageLength(dist);

        // Checks message is correct
        auto const sendLengths = xc.GetSendLengths();
        auto const neighbors = graph.GetNeighbors();
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), sendLengths.size());
        for(auto const item: util::zip(neighbors, sendLengths))
        {
          auto const sending = std::get<0>(item) == sendto;
          size_t const nCells = sending ? 1: 0;
          size_t const nVertices = sending ? (*cells.begin())->GetNumberOfNodes(): 0;
          CPPUNIT_ASSERT_EQUAL(nCells, std::get<1>(item).nCells);
          CPPUNIT_ASSERT_EQUAL(nVertices, std::get<1>(item).nVertices);
        }

        // Wait for end of request and check received lengths
        xc.WaitOutLengthRequest();
        auto const recvfrom =
          graph.Rank() == 0 ? 1:
          graph.Rank() == 1 ? 2:
          graph.Rank() == 2 ? 3: std::numeric_limits<size_t>::max();

        auto const recvLengths = xc.GetReceiveLengths();
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), recvLengths.size());
        for(auto const item: util::zip(neighbors, recvLengths))
        {
          auto const receiving = std::get<0>(item) == recvfrom;
          size_t const nCells = receiving ? 1: 0;
          size_t const nVerts = receiving ? GetCell(center, 1e0, recvfrom)->GetNumberOfNodes(): 0;
          CPPUNIT_ASSERT_EQUAL(nCells, std::get<1>(item).nCells);
          CPPUNIT_ASSERT_EQUAL(nVerts, std::get<1>(item).nVertices);
        }
      }

      CPPUNIT_TEST_SUITE_REGISTRATION (CellParallelizationTests);
    }
  }
}

#endif  // ONCE
