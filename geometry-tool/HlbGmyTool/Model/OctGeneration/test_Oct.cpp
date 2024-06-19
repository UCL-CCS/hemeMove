// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <iostream>
#include <memory>

#include "Oct.h"

namespace hemelb::gmytool::oct {

struct OctreeTests {
  using Tree = Octree<char>;
  using TreePtr = std::shared_ptr<Tree>;

  TreePtr mk_8cube_234() {
    auto tree = std::make_shared<Tree>(3);
    auto node = tree->GetCreate(2, 3, 4, 0);
    node->Data() = 1;
    return tree;
  }

  void CreateEmpty() { auto tree = Tree(3); }

  void CreateSimple() {
    auto tree = mk_8cube_234();
    for (auto i : {0, 4})
      for (auto j : {0, 4})
        for (auto k : {0, 4}) {
          // Most of these should be nullptr
          auto tmp = tree->Get(i, j, k, 2);
          if (i == 0 && j == 0 && k == 4)
            // But not this one
            REQUIRE(tmp);
          else
            REQUIRE(!tmp);
        }

    auto n2 = tree->Get(2, 2, 4, 1);
    REQUIRE(n2);
    auto n3 = tree->Get(2, 3, 4, 0);
    REQUIRE(n3);
    REQUIRE(n3->Data() == 1);
  }

  class Counter : public Tree::Visitor {
   public:
    int i;

    Counter() : i(0) {}

    virtual void Do(Tree::Node& node) {
      REQUIRE(node.Level() == i);
      ++i;
    }
  };

  void SimpleIter() {
    auto tree = mk_8cube_234();
    Tree::Int expected_nodes[4][3] = {
        {2, 3, 4}, {2, 2, 4}, {0, 0, 4}, {0, 0, 0}};
    // iterate down through all nodes
    int i = 0;
    tree->IterDepthFirst([&i, &expected_nodes](Tree::NodePtr node) mutable {
      REQUIRE(node->Level() == i);
      REQUIRE(node->X() == expected_nodes[i][0]);
      REQUIRE(node->Y() == expected_nodes[i][1]);
      REQUIRE(node->Z() == expected_nodes[i][2]);
      ++i;
    });
    REQUIRE(i == 4);

    i = 1;
    tree->IterDepthFirst(1, [&i](Tree::NodePtr node) mutable {
      REQUIRE(node->Level() == i);
      ++i;
    });
    REQUIRE(i == 4);
  }
};

METHOD_AS_TEST_CASE(OctreeTests::CreateEmpty, "CreateEmpty", "[oct]");
METHOD_AS_TEST_CASE(OctreeTests::CreateSimple, "CreateSimple", "[oct]");
METHOD_AS_TEST_CASE(OctreeTests::SimpleIter, "SimpleIter", "[oct]");

}  // namespace hemelb::gmytool::oct
