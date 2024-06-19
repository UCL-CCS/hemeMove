// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <set>

#include "TestResources/Meshes.h"
#include "TriangleSorter.h"

namespace hemelb::gmytool::oct {

struct TriangleSorterTests {
  void Trivial() {
    // 16 cube
    TriTree::Int levels = 4;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();
    // Sort!
    auto tree =
        TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);

    REQUIRE(tree.Level() == levels);
    TriTree::Int i = tri_level;
    tree.IterDepthFirst([&i](TriTree::NodePtr node) mutable {
      REQUIRE(node->X() == 0);
      REQUIRE(node->Y() == 0);
      REQUIRE(node->Z() == 0);
      REQUIRE(node->Level() == i);
      ++i;
    });

    auto node = tree.Get(0, 0, 0, tri_level);
    auto ids = node->Data();
    std::sort(ids.begin(), ids.end());
    REQUIRE(ids.size() == 2);
    REQUIRE(ids.find(0) != ids.end());
    REQUIRE(ids.find(1) != ids.end());
  }

  void TrivialDeep() {
    // 1024 cube
    TriTree::Int levels = 10;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();
    // Sort!
    auto tree =
        TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);

    REQUIRE(tree.Level() == levels);
    TriTree::Int i = tri_level;
    tree.IterDepthFirst([&i](TriTree::NodePtr node) mutable {
      REQUIRE(node->X() == 0);
      REQUIRE(node->Y() == 0);
      REQUIRE(node->Z() == 0);
      REQUIRE(node->Level() == i);
      ++i;
    });

    auto node = tree.Get(0, 0, 0, tri_level);
    auto ids = node->Data();
    std::sort(ids.begin(), ids.end());
    REQUIRE(ids.size() == 2);
    REQUIRE(ids.find(0) != ids.end());
    REQUIRE(ids.find(1) != ids.end());
  }

  template <class T>
  bool overlap1d(T aMin, T aMax, T bMin, T bMax) {
    return (aMin < bMax) && (aMax > bMin);
  }

  bool overlap3d(const Vector& aMin,
                 const Vector& aMax,
                 const Vector& bMin,
                 const Vector& bMax) {
    bool ans = true;
    for (auto d = 0; d < 3; ++d)
      ans &= overlap1d(aMin[d], aMax[d], bMin[d], bMax[d]);

    return ans;
  }

  void Sphere() {
    TriTree::Int levels = 5;
    TriTree::Int tri_level = 3;

    auto sphere = SimpleMeshFactory::MkSphere();
    auto tree = TrianglesToTreeSerial(levels, tri_level, sphere->points,
                                      sphere->triangles);

    std::set<int> seen_tris;

    tree.IterDepthFirst(tri_level, tri_level, [&](TriTree::NodePtr node) {
      const auto& triIds = node->Data();
      const Vector offset(node->X(), node->Y(), node->Z());
      auto node_min = offset - Vector{1};
      auto node_max = offset + Vector{(1 << tri_level) + 1};

      for (auto iTri : triIds) {
        Vector tri_min(std::numeric_limits<double>::infinity());
        Vector tri_max(0);
        for (auto iPt : sphere->triangles[iTri]) {
          tri_min.UpdatePointwiseMin(sphere->points[iPt]);
          tri_max.UpdatePointwiseMax(sphere->points[iPt]);
        }

        REQUIRE(overlap3d(tri_min, tri_max, node_min, node_max));
        seen_tris.insert(triIds.begin(), triIds.end());
      }
    });

    // Assert that we have seen all the triangles in the input
    REQUIRE(seen_tris.size() == sphere->triangles.size());
    int i = 0;
    for (auto seen_id : seen_tris) {
      REQUIRE(seen_id == i);
      ++i;
    }
  }

  bool trees_with_triIds_equal(TriTree& t1,
                               TriTree& t2,
                               TriTree::Int tri_level) {
    bool eq = true;
    try {
      t1.IterDepthFirst(tri_level, tri_level, [&](TriTree::NodePtr n1) {
        auto n2 = t2.Get(n1->X(), n1->Y(), n1->Z(), n1->Level());
        if (n1->Data() != n2->Data())
          throw 0;
      });
    } catch (int& i) {
      eq = false;
    }
    return eq;
  }
  void MergeIdentity() {
    // 16 cube
    TriTree::Int levels = 4;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();

    auto tree =
        TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    TriTree nulltree(levels);

    TreeSummer s(levels, tri_level);
    s.Add(tree);
    s.Add(nulltree);
    auto newtree = s.GetTree();
    REQUIRE(trees_with_triIds_equal(newtree, tree, tri_level));
  }

  void MergeSimple() {
    // 16
    TriTree::Int levels = 4;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();

    auto start = triv->triangles.cbegin();
    auto middle = start + 1;
    auto end = triv->triangles.cend();

    auto t1 = TrianglesToTree_Worker(levels, tri_level, triv->points, start,
                                     middle, 0);
    auto t2 =
        TrianglesToTree_Worker(levels, tri_level, triv->points, middle, end, 1);

    TreeSummer s(levels, tri_level);
    s.Add(t1);
    s.Add(t2);
    auto newtree = s.GetTree();

    auto tree =
        TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    REQUIRE(trees_with_triIds_equal(newtree, tree, tri_level));
  }

  void ParallelTrivial() {
    // 16 cube
    TriTree::Int levels = 4;
    // put triangles onto the 8 cube level
    TriTree::Int tri_level = 3;
    auto triv = SimpleMeshFactory::MkTrivial();
    // Sort!
    auto t1 =
        TrianglesToTreeSerial(levels, tri_level, triv->points, triv->triangles);
    auto t2 = TrianglesToTreeParallel(levels, tri_level, triv->points,
                                      triv->triangles, 2);
    REQUIRE(trees_with_triIds_equal(t1, t2, tri_level));
  }

  void ParallelSphere() {
    int np = GENERATE(0, 2, 3, 4);
    TriTree::Int levels = 5;
    TriTree::Int tri_level = 3;

    auto sphere = SimpleMeshFactory::MkSphere();
    auto t1 = TrianglesToTreeSerial(levels, tri_level, sphere->points,
                                    sphere->triangles);
    auto t2 = TrianglesToTreeParallel(levels, tri_level, sphere->points,
                                      sphere->triangles, np);
    REQUIRE(trees_with_triIds_equal(t1, t2, tri_level));
  }

  void Duct() {
    // 16 cube
    TriTree::Int levels = 4;
    TriTree::Int tri_level = 2;
    auto duct = SimpleMeshFactory::MkDuct();
    auto tree =
        TrianglesToTreeSerial(levels, tri_level, duct->points, duct->triangles);

    REQUIRE(tree.Get(0, 0, 0, 3));
    REQUIRE(tree.Get(0, 0, 8, 3));

    REQUIRE(tree.Get(0, 0, 0, 2));
    REQUIRE(tree.Get(0, 0, 4, 2));
    REQUIRE(tree.Get(0, 4, 0, 2));
    REQUIRE(tree.Get(0, 4, 4, 2));
    REQUIRE(tree.Get(4, 0, 0, 2));
    REQUIRE(tree.Get(4, 0, 4, 2));
    REQUIRE(tree.Get(4, 4, 0, 2));
    REQUIRE(tree.Get(4, 4, 4, 2));
  }
};

METHOD_AS_TEST_CASE(TriangleSorterTests::Trivial, "Trivial", "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::TrivialDeep,
                    "TrivialDeep",
                    "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::Sphere, "Sphere", "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::MergeIdentity,
                    "MergeIdentity",
                    "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::MergeSimple,
                    "MergeSimple",
                    "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::ParallelTrivial,
                    "ParallelTrivial",
                    "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::ParallelSphere,
                    "ParallelSphere",
                    "[triangle]");
METHOD_AS_TEST_CASE(TriangleSorterTests::Duct, "Duct", "[triangle]");

}  // namespace hemelb::gmytool::oct
