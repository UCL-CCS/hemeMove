//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_NET_NEIGHBORSCOMMTESTS_H
#define HEMELB_UNITTESTS_NET_NEIGHBORSCOMMTESTS_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <functional>

#include "net/MpiCommunicator.h"
#include "net/INeighborAllToAll.h"
#include "net/INeighborAllToAllV.h"
#include "util/Iterator.h"

namespace hemelb
{
  namespace unittests
  {
    namespace net
    {
      using namespace hemelb::net;
      class NeighborCommTests : public CppUnit::TestFixture
      {
          CPPUNIT_TEST_SUITE (NeighborCommTests);
          CPPUNIT_TEST(testInteger);
          CPPUNIT_TEST(testDoubles);
          CPPUNIT_TEST(testVariableInts);
          CPPUNIT_TEST_SUITE_END();

        public:
          void setUp();
          //! Sends one integer per proc
          void testInteger();
          //! Sends several doubles per proc
          void testDoubles();
          //! \brief Sends different number of ints from each proc
          //! \details Each proc expects to know how many it is gettin from each before hand
          void testVariableInts();

        protected:
          net::MpiCommunicator graph;
      };

      void NeighborCommTests::setUp()
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

      void NeighborCommTests::testInteger()
      {
        if(not graph)
        {
          return;
        }
        // Just a message such what i sends to j is different from what j sends to i
        auto const message = [](int proci, int procj)
        {
          return proci + procj + (proci > procj ? 10: proci < procj ? -10: 0);
        };

        // Prepares buffer. One message per object.
        auto const neighbors = graph.GetNeighbors();
        INeighborAllToAll<int> all2all(graph);
        all2all.GetSendBuffer().resize(neighbors.size());
        for(auto const item: util::zip(all2all.GetSendBuffer(), neighbors))
        {
          std::get<0>(item) = message(graph.Rank(), std::get<1>(item));
        }
        // Send message
        all2all.send();
        // Wait for end of request and check received lengths
        all2all.receive();

        CPPUNIT_ASSERT_EQUAL(neighbors.size(), all2all.GetReceiveBuffer().size());
        for(auto const item: util::zip(all2all.GetReceiveBuffer(), neighbors))
        {
          auto const received = std::get<0>(item);
          auto const neighbor = std::get<1>(item);
          CPPUNIT_ASSERT_EQUAL(message(neighbor, graph.Rank()), received);
        }
      }

      void NeighborCommTests::testDoubles()
      {
        if(not graph)
        {
          return;
        }
        // Just a message such what i sends to j is different from what j sends to i
        auto const message = [](int proci, int procj, int n = 0)
        {
          return n + proci + procj + (proci > procj ? 100: proci < procj ? -100: 0);
        };

        // Prepares buffer. One message per object.
        auto const neighbors = graph.GetNeighbors();
        INeighborAllToAll<int> all2all(graph);
        for(auto const n: neighbors)
        {
          auto const r = graph.Rank();
          all2all.AddToBuffer(n, {message(r, n), message(r, n, 1), message(r, n, 2)});
        }
        // Send message
        all2all.send();
        // Wait for end of request and check received lengths
        all2all.receive();

        auto const & received = all2all.GetReceiveBuffer();
        CPPUNIT_ASSERT_EQUAL(neighbors.size() * 3, received.size());
        for(auto const item: util::enumerate(neighbors))
        {
          CPPUNIT_ASSERT_EQUAL(message(item.value, graph.Rank()), received[item.index * 3]);
          CPPUNIT_ASSERT_EQUAL(message(item.value, graph.Rank(), 1), received[item.index * 3 + 1]);
          CPPUNIT_ASSERT_EQUAL(message(item.value, graph.Rank(), 2), received[item.index * 3 + 2]);
        }
      }

      void NeighborCommTests::testVariableInts()
      {
        if(not graph)
        {
          return;
        }
        // Just a message such what i sends to j is different from what j sends to i
        auto const message = [](int proci, int procj, int n = 0)
        {
          return n + proci + procj + (proci > procj ? 100: proci < procj ? -100: 0);
        };
        // If X sends to Y, then how many it does send
        auto const XSendsToY = [](int proci, int procj)
        {
          return proci >= 4 or procj >= 4 ? 0:
              proci > procj ? proci + procj: proci + procj + 10;
        };

        // Prepares arrays how many elements are sent from each to each
        auto const neighbors = graph.GetNeighbors();
        std::vector<int> sendCounts, recvCounts;
        int nReceived = 0;
        for(auto neighbor: neighbors)
        {
          sendCounts.push_back(XSendsToY(graph.Rank(), neighbor));
          recvCounts.push_back(XSendsToY(neighbor, graph.Rank()));
          nReceived += recvCounts.back();
        }

        // Prepares array with elements to send
        INeighborAllToAllV<int> all2all(graph, sendCounts, recvCounts);
        for(auto const neighbor: neighbors)
        {
          std::vector<int> stuff;
          for(int i(0); i < XSendsToY(graph.Rank(), neighbor); ++i)
          {
            stuff.push_back(message(graph.Rank(), neighbor, i));
          }
          all2all.insertSend(neighbor, stuff);
        }

        // Send message
        all2all.send();
        // Wait for end of request and check received lengths
        all2all.receive();

        // Check result
        CPPUNIT_ASSERT_EQUAL(size_t(nReceived), all2all.GetReceiveBuffer().size());
        auto received = all2all.GetReceiveBuffer().begin();
        for(auto const item: util::zip(neighbors, recvCounts))
        {
          for(int i(0); i < std::get<1>(item); ++i, ++received)
          {
            CPPUNIT_ASSERT_EQUAL(message(std::get<0>(item), graph.Rank(), i), *received);
          }
        }

      }

      CPPUNIT_TEST_SUITE_REGISTRATION (NeighborCommTests);
    }
  }
}

#endif  // ONCE

