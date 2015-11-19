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

#         define HEMELB_MACRO(Name, name, TYPE)   \
            net::TYPE & Get ## Name()             \
            {                                     \
              return name;                        \
            }

            HEMELB_MACRO(CellCount, cellCount, INeighborAllToAll<int>);
            HEMELB_MACRO(TotalNodeCount, totalNodeCount, INeighborAllToAll<int>);
            HEMELB_MACRO(NameLengths, nameLengths, INeighborAllToAll<int>);
            HEMELB_MACRO(NodeCount, nodeCount, INeighborAllToAllV<size_t>);
            HEMELB_MACRO(CellScales, cellScales, INeighborAllToAllV<LatticeDistance>);
            HEMELB_MACRO(OwnerIDs, ownerIDs, INeighborAllToAllV<int>);
            HEMELB_MACRO(CellUUIDs, cellUUIDs, INeighborAllToAllV<unsigned char>);
            HEMELB_MACRO(NodePositions, nodePositions, INeighborAllToAllV<LatticeDistance>);
            HEMELB_MACRO(TemplateNames, templateNames, INeighborAllToAllV<char>);
#         undef HEMELB_MACRO
      };

      using namespace hemelb::redblood;
      using namespace hemelb::redblood::parallel;
      class CellParallelizationTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (CellParallelizationTests);
          CPPUNIT_TEST(testCellSwapGetLength);
          CPPUNIT_TEST(testCellSwapPostCells);
          CPPUNIT_TEST(testSingleCellSwapWithRetainedOwnership);
          CPPUNIT_TEST(testSingleCellSwap);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();
          //! Test message sending number of cells and number of nodes
          void testCellSwapGetLength();
          //! Test messages from swapping cells
          void testCellSwapPostCells();
          //! Test messages from swapping cells while retaining ownership
          void testSingleCellSwapWithRetainedOwnership();
          //! Test messages from swapping cells while retaining ownership
          void testSingleCellSwap();

          //! Set of nodes affected by given proc
          std::set<proc_t> nodeLocation(LatticePosition const &node);

          //! Center of the region affected by a node
          LatticePosition GetCenter(proc_t rank) const;

          //! Creates map of distributed nodes
          CellParallelization::NodeDistributions GetNodeDistribution(
              CellContainer const &cells) const;

          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(
              LatticePosition const &pos, boost::uuids::uuid const & uuid,
              Dimensionless scale=1e0, unsigned int depth=0) const;
          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(
              LatticePosition const &pos, Dimensionless scale=1e0, unsigned int depth=0) const
          {
            return GetCell(pos, boost::uuids::uuid(), scale, depth);
          }
          std::shared_ptr<Cell> GivenCell(size_t i) const;
          //! Get cell centered at given position
          std::shared_ptr<Cell> GetCell(Dimensionless scale=1e0, unsigned int depth=0) const
          {
            return GetCell(GetCenter(graph ? graph.Rank(): 0), scale, depth);
          }

          //! Id of owning cell
          int Ownership(CellContainer::const_reference cell) const;

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
          LatticePosition const &pos, boost::uuids::uuid const & uuid,
          Dimensionless scale, unsigned int depth) const
      {
        auto cell = std::make_shared<Cell>(icoSphere(depth));
        *cell += pos - cell->GetBarycenter();
        *cell *= scale;
        cell->SetScale(scale);
        cell->SetTag(uuid);
        return cell;
      }

      std::shared_ptr<Cell> CellParallelizationTests::GivenCell(size_t process) const
      {
        size_t const sendto =
          process > 0 and process < 4 ? process - 1: std::numeric_limits<size_t>::max();
        auto const center = GetCenter(sendto);
        auto const scale = 1.0 + 0.1 * static_cast<double>(process);
        boost::uuids::uuid uuid;
        std::fill(uuid.begin(), uuid.end(), static_cast<char>(process));
        auto result = GetCell(center, uuid, scale, process);
        std::ostringstream sstr; sstr << process;
        result->SetTemplateName(sstr.str());
        return result;
      };

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

      int CellParallelizationTests::Ownership(CellContainer::const_reference cell) const
      {
        std::vector<LatticeDistance> distances;
        auto const ranks = graph.RankMap(net::MpiCommunicator::World());
        auto const barycenter = cell->GetBarycenter();
        for(size_t i(0); i < graph.Size(); ++i)
        {
          auto const d = (GetCenter(ranks.find(i)->second) - barycenter).GetMagnitudeSquared();
          distances.push_back(d);
        };
        return std::min_element(distances.begin(), distances.end()) - distances.begin();
      }

      // Check that cell send each other whole cells
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
        CPPUNIT_ASSERT(xc.GetCellCount().GetCommunicator());
        CPPUNIT_ASSERT(xc.GetTotalNodeCount().GetCommunicator());
        xc.PostCellMessageLength(dist, cells);

        // Checks message is correct
        auto const &sendCellCount = xc.GetCellCount().GetSendBuffer();
        auto const &sendTotalNodeCount = xc.GetTotalNodeCount().GetSendBuffer();
        auto const neighbors = graph.GetNeighbors();
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), sendCellCount.size());
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), sendTotalNodeCount.size());
        for(auto const item: util::zip(neighbors, sendCellCount, sendTotalNodeCount))
        {
          auto const sending = std::get<0>(item) == sendto;
          size_t const nCells = sending ? 1: 0;
          size_t const nVertices = sending ? (*cells.begin())->GetNumberOfNodes(): 0;
          CPPUNIT_ASSERT_EQUAL(int(nCells), std::get<1>(item));
          CPPUNIT_ASSERT_EQUAL(int(nVertices), int(std::get<2>(item)));
        }

        // Wait for end of request and check received lengths
        xc.GetCellCount().receive();
        xc.GetTotalNodeCount().receive();
        xc.GetNameLengths().receive();
        auto const recvfrom =
          graph.Rank() == 0 ? 1:
          graph.Rank() == 1 ? 2:
          graph.Rank() == 2 ? 3: std::numeric_limits<size_t>::max();

        auto const &receiveCellCount = xc.GetCellCount().GetReceiveBuffer();
        auto const &receiveTotalNodeCount = xc.GetTotalNodeCount().GetReceiveBuffer();
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), receiveCellCount.size());
        CPPUNIT_ASSERT_EQUAL(neighbors.size(), receiveTotalNodeCount.size());
        for(auto const item: util::zip(neighbors, receiveCellCount, receiveTotalNodeCount))
        {
          auto const receiving = std::get<0>(item) == recvfrom;
          size_t const nCells = receiving ? 1: 0;
          size_t const nVerts = receiving ? GetCell(center, 1e0, recvfrom)->GetNumberOfNodes(): 0;
          CPPUNIT_ASSERT_EQUAL(int(nCells), std::get<1>(item));
          CPPUNIT_ASSERT_EQUAL(int(nVerts), int(std::get<2>(item)));
        }
      }

      void CellParallelizationTests::testCellSwapPostCells()
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
        auto const getScale = [](int process)
        {
          return 1.0 + 0.1 * static_cast<double>(process);
        };
        CellContainer cells{GetCell(center, getScale(graph.Rank()), graph.Rank())};
        if(graph.Rank() == 0)
        {
          cells.clear();
        }
        auto const dist = GetNodeDistribution(cells);

        ExchangeCells xc(graph);
        xc.PostCellMessageLength(dist, cells);

        auto keepOwnership = [](CellContainer::const_reference)
        {
          return net::MpiCommunicator::World();
        };
        xc.PostCells(dist, cells, keepOwnership);

        // check message sizes
        auto const neighbors = graph.GetNeighbors();
        unsigned long const Nsend = graph.Rank() > 0 and graph.Rank() < 4 ? 1: 0;
        unsigned long const Nreceive = graph.Rank() < 3 ? 1: 0;
        unsigned long const uuid_size = sizeof(boost::uuids::uuid);
        CPPUNIT_ASSERT_EQUAL(Nsend, xc.GetCellScales().GetSendBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nreceive, xc.GetCellScales().GetReceiveBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nsend, xc.GetNodeCount().GetSendBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nreceive, xc.GetNodeCount().GetReceiveBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nsend * uuid_size , xc.GetCellUUIDs().GetSendBuffer().size());
        CPPUNIT_ASSERT_EQUAL(Nreceive * uuid_size , xc.GetCellUUIDs().GetReceiveBuffer().size());

        if(cells.size() > 0)
        {
          auto const cell = *cells.begin();
          auto const xcScale = xc.GetCellScales().GetSendBuffer()[0];
          auto const xcTag = xc.GetCellUUIDs().GetSendBuffer()[0];
          auto const xcNodes = xc.GetNodeCount().GetSendBuffer()[0];
          CPPUNIT_ASSERT_DOUBLES_EQUAL(cell->GetScale(), xcScale, 1e-8);
          CPPUNIT_ASSERT_EQUAL(*cell->GetTag().begin(), xcTag);
          CPPUNIT_ASSERT_EQUAL(cell->GetNumberOfNodes(), site_t(xcNodes));
        }

        // receive messages
        xc.GetNodeCount().receive();
        xc.GetCellScales().receive();
        xc.GetOwnerIDs().receive();
        xc.GetCellUUIDs().receive();
        xc.GetNodePositions().receive();
        xc.GetTemplateNames().receive();

        if(graph.Rank() < 3)
        {
          auto const scale = getScale(graph.Rank() + 1);
          auto const nNodes = GetCell(center, scale, graph.Rank() + 1)->GetNumberOfNodes();
          CPPUNIT_ASSERT_DOUBLES_EQUAL(scale, xc.GetCellScales().GetReceiveBuffer()[0], 1e-8);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(nNodes, xc.GetNodeCount().GetReceiveBuffer()[0], 1e-8);
        }
      }

      void CellParallelizationTests::testSingleCellSwapWithRetainedOwnership()
      {
        if(not graph)
        {
          return;
        }

        auto templates = std::make_shared<TemplateCellContainer>(
            TemplateCellContainer{{"1", GivenCell(1)}, {"2", GivenCell(2)}, {"3", GivenCell(3)}});

        CellContainer owned{GivenCell(graph.Rank())};
        std::map<proc_t, CellContainer> lent;
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        xc.PostCellMessageLength(dist, owned);
        auto keepOwnership = [](CellContainer::const_reference)
        {
          return net::MpiCommunicator::World().Rank();
        };
        xc.PostCells(dist, owned, keepOwnership);
        xc.ReceiveCells(owned, lent, templates);

        if(graph.Rank() < 3)
        {
          CPPUNIT_ASSERT_EQUAL(size_t(1), lent.size());
          CPPUNIT_ASSERT_EQUAL(size_t(1), lent.begin()->second.size());
          HEMELB_CAPTURE2(proc_t(graph.Rank() + 1), lent.begin()->first);
          CPPUNIT_ASSERT_EQUAL(proc_t(graph.Rank() + 1), lent.begin()->first);
          auto const actual = *lent.begin()->second.begin();
          auto const expected = GivenCell(graph.Rank()+1);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(expected->GetScale(), actual->GetScale(), 1e-8);
          CPPUNIT_ASSERT_EQUAL(expected->GetTemplateName(), actual->GetTemplateName());
          CPPUNIT_ASSERT(expected->GetTag() == actual->GetTag());
          CPPUNIT_ASSERT_EQUAL(expected->GetNumberOfNodes(), actual->GetNumberOfNodes());
          for(auto item: util::zip(expected->GetVertices(), actual->GetVertices()))
          {
            CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
          }
        }
        else
        {
          CPPUNIT_ASSERT_EQUAL(size_t(0), lent.size());
        }
      }

      void CellParallelizationTests::testSingleCellSwap()
      {
        if(not graph)
        {
          return;
        }

        auto templates = std::make_shared<TemplateCellContainer>(
            TemplateCellContainer{{"1", GivenCell(1)}, {"2", GivenCell(2)}, {"3", GivenCell(3)}});

        CellContainer owned{GivenCell(graph.Rank())};
        std::map<proc_t, CellContainer> lent;
        auto const dist = GetNodeDistribution(owned);

        ExchangeCells xc(graph);
        xc.PostCellMessageLength(dist, owned);
        auto const ownership = [this](CellContainer::const_reference cell)
        {
          return Ownership(cell);
        };
        xc.PostCells(dist, owned, ownership);
        xc.ReceiveCells(owned, lent, templates);

        CPPUNIT_ASSERT_EQUAL(size_t(0), lent.size());
        // if(graph.Rank() < 3)
        // {
        //   CPPUNIT_ASSERT_EQUAL(size_t(1), lent.size());
        //   CPPUNIT_ASSERT_EQUAL(size_t(1), lent.begin()->second.size());
        //   HEMELB_CAPTURE2(proc_t(graph.Rank() + 1), lent.begin()->first);
        //   CPPUNIT_ASSERT_EQUAL(proc_t(graph.Rank() + 1), lent.begin()->first);
        //   auto const actual = *lent.begin()->second.begin();
        //   auto const expected = GivenCell(graph.Rank()+1);
        //   CPPUNIT_ASSERT_DOUBLES_EQUAL(expected->GetScale(), actual->GetScale(), 1e-8);
        //   CPPUNIT_ASSERT_EQUAL(expected->GetTemplateName(), actual->GetTemplateName());
        //   CPPUNIT_ASSERT(expected->GetTag() == actual->GetTag());
        //   CPPUNIT_ASSERT_EQUAL(expected->GetNumberOfNodes(), actual->GetNumberOfNodes());
        //   for(auto item: util::zip(expected->GetVertices(), actual->GetVertices()))
        //   {
        //     CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).x, std::get<1>(item).x, 1e-8);
        //     CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).y, std::get<1>(item).y, 1e-8);
        //     CPPUNIT_ASSERT_DOUBLES_EQUAL(std::get<0>(item).z, std::get<1>(item).z, 1e-8);
        //   }
        // }
        // else
        // {
        // }
      }
      CPPUNIT_TEST_SUITE_REGISTRATION (CellParallelizationTests);
    }
  }
}

#endif  // ONCE
