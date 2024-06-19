// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <map>
#include <set>

#include "H5.h"
#include "enumerate.hpp"

namespace hemelb::gmytool::oct {

// Check our iterator types satisfy the concept.
//
// Per DR P2325R3, std::weakly_incrementable (required by iterator
// concept) should NOT require default constructible, but not all std
// library implementations have caught up :(
namespace {
struct NonDefaultConstructible {
  using difference_type = int;
  using value_type = int;
  using reference = int&;
  using pointer = int*;
  using iterator_category = std::input_iterator_tag;

  NonDefaultConstructible(int);
  NonDefaultConstructible& operator++();
  NonDefaultConstructible operator++(int);
};

constexpr bool DEFAULT_CONSTRUCTIBLE_REQUIRED =
    !std::weakly_incrementable<NonDefaultConstructible>;
}  // namespace
static_assert(DEFAULT_CONSTRUCTIBLE_REQUIRED ||
              std::input_iterator<H5::LinkIterator>);
static_assert(DEFAULT_CONSTRUCTIBLE_REQUIRED ||
              std::input_iterator<H5::AttrIterator>);

auto mk_empty_in_memory_file() {
  auto fapl = H5::PList::FileAccess();
  H5Pset_fapl_core(*fapl, 1 << 20, false);
  return H5::File::Create("inmem.h5", H5F_ACC_TRUNC, H5::PList::Default(),
                          fapl);
}

auto const some_names = std::vector<std::string>{"a", "b", "c", "d", "e"};

TEST_CASE("Basic group manipulation", "[h5]") {
  auto file = mk_empty_in_memory_file();
  auto foo = file->CreateGroup("foo");

  for (auto&& n : some_names) {
    foo->CreateGroup(n);
  }

  SECTION("Can iterate groups") {
    auto unseen_names =
        std::set<std::string>{some_names.begin(), some_names.end()};
    int i = 0;

    for (auto&& g_name : *foo) {
      REQUIRE(unseen_names.contains(g_name));
      unseen_names.erase(g_name);
    }
    REQUIRE(unseen_names.size() == 0);
  }
}

TEST_CASE("Basic attributes work", "[h5]") {
  auto file = mk_empty_in_memory_file();
  auto root = file->CreateGroup("root");
  int i = 0;
  std::map<std::string, int> expected_attrs;
  for (auto&& n : some_names) {
    root->SetAttribute(n, i);
    expected_attrs[n] = i;
    ++i;
  }
  REQUIRE(i == some_names.size());

  auto iterhelp = H5::Attribute::Iterate(*root);
  for (auto&& key : iterhelp) {
    REQUIRE(expected_attrs.contains(key));
    int val;
    root->GetAttribute(key, val);
    REQUIRE(expected_attrs[key] == val);
    expected_attrs.erase(key);
  }
  REQUIRE(expected_attrs.size() == 0);
}

TEST_CASE("More complex attributes work", "[h5]") {
  auto file = mk_empty_in_memory_file();
  auto root = file->CreateGroup("root");
  root->SetAttribute("int", 101);
  root->SetAttribute("float", -4.f);
  std::vector<double> data = {0.0, 1.0, 2.0, 3.0};
  root->SetAttribute("vector", data);
  std::string text = "This is some textual data to be stored";
  root->SetAttribute("string", text);

  // Now pull them out
  {
    int i;
    root->GetAttribute("int", i);
    REQUIRE(i == 101);
  }
  {
    float f;
    root->GetAttribute("float", f);
    REQUIRE(f == -4.f);
  }
  {
    std::vector<double> actual;
    root->GetAttribute("vector", actual);
    for (auto [i, x] : enumerate(data)) {
      REQUIRE(x == actual[i]);
    }
  }
  {
    std::string actual;
    root->GetAttribute("string", actual);
    REQUIRE(actual == text);
  }
  {
    int i;
    REQUIRE_THROWS_AS(root->GetAttribute("Nonexisting", i), H5::Error);
    REQUIRE_THROWS_AS(root->GetAttribute("string", i), H5::Error);
  }
}

}  // namespace hemelb::gmytool::oct
