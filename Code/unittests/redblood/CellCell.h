//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLCELL_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLCELL_H

#include <cppunit/TestFixture.h>
#include "redblood/Cell.h"
#include "unittests/redblood/Fixtures.h"
#include "unittests/redblood/CellCell.h"
#include "redblood/CellCell.cc"

namespace hemelb { namespace unittests { namespace redblood {

class CellCellInteractionTests : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(CellCellInteractionTests);
     CPPUNIT_TEST(testBoxHaloTooBig);
     CPPUNIT_TEST(testBoxHalo);
     CPPUNIT_TEST(testAddNodes);
     CPPUNIT_TEST(testAddMeshes);
     CPPUNIT_TEST(testIterator);
     CPPUNIT_TEST(testUpdate);
    CPPUNIT_TEST_SUITE_END();

public:
    void testBoxHaloTooBig();
    void testBoxHalo();
    void testAddNodes();
    void testAddMeshes();
    void testIterator();
    void testUpdate();

protected:
    boost::shared_ptr<CellContainer> fixture(PhysicalDistance) const;
};

void CellCellInteractionTests :: testBoxHaloTooBig() {
  PhysicalDistance const cutoff(5.0);
  DivideConquer<int> dnc(cutoff);
  LatticeVector const key(0, 1, -2);
  LatticePosition const position
    = (LatticePosition(key) + LatticePosition(0.6, 0.5, 0.5)) * cutoff;

  CPPUNIT_ASSERT(not figure_nearness(dnc, key, position, cutoff * 0.501));
  CPPUNIT_ASSERT(figure_nearness(dnc, key, position, cutoff * 0.499));
}

void CellCellInteractionTests :: testBoxHalo() {
  DivideConquer<int> dnc(5.0);
  LatticeVector const key(0, 1, -2);
  LatticePosition const center =
    (LatticePosition(key) + LatticePosition(0.5, 0.5, 0.5)) * dnc.GetBoxSize();
  CPPUNIT_ASSERT(figure_nearness(dnc, key, center, 2.0) == 0);
  for(size_t d(1); d < (1 << 6); d <<= 1) {
    LatticePosition const disp = CellReference::directions(d) * 0.6;
    int const nearness = figure_nearness(dnc, key, center + disp, 2.0);
    CPPUNIT_ASSERT(nearness == d);
  }

  LatticePosition const mult =
    CellReference::directions(CellReference::TOP) * 0.6
    + CellReference::directions(CellReference::NORTH) * 0.6
    + CellReference::directions(CellReference::EAST) * 0.6;
  int const actual = figure_nearness(dnc, key, center + mult, 2.0);
  int const expected
    = CellReference::TOP | CellReference::NORTH | CellReference::EAST;
  CPPUNIT_ASSERT(actual == expected);
}

void CellCellInteractionTests :: testAddNodes() {
  PhysicalDistance const cutoff = 5.0;
  PhysicalDistance const halo = 2.0;
  site_t const cellIndex(4);
  DivideConquer<CellReference> dnc(cutoff);

  // Adds nodes, last in halo
  LatticePosition const center = LatticePosition(1, 1, 1) * (0.5 * cutoff);
  PhysicalDistance const offhalo = (cutoff * 0.5 - 1.1 * halo) / cutoff;
  PhysicalDistance const inhalo = (cutoff * 0.5 - 0.9 * halo) / cutoff;
  MeshData::t_Vertices vertices;
  vertices.push_back(center);
  vertices.push_back(center + LatticePosition(offhalo, 0, 0) * cutoff);
  vertices.push_back(center + LatticePosition(2, offhalo + 3.0, -2) * cutoff);
  vertices.push_back(
      center + LatticePosition(1, 0, 0) * cutoff
      + CellReference::directions(CellReference::NORTH) * inhalo * cutoff
  );

  initialize_cells(dnc, vertices, cellIndex, halo);
  CPPUNIT_ASSERT(dnc.size() == vertices.size());

  DivideConquer<CellReference>::const_range const omega
    = dnc.equal_range(LatticeVector(0, 0, 0));
  DivideConquer<CellReference>::const_range const empty
    = dnc.equal_range(LatticeVector(10, 10, 10));
  DivideConquer<CellReference>::const_range const alpha
    = dnc.equal_range(LatticeVector(2, 3, -2));
  DivideConquer<CellReference>::const_range const haloed
    = dnc.equal_range(LatticeVector(1, 0, 0));

  CPPUNIT_ASSERT(std::distance(omega.first, omega.second) == 2);
  CPPUNIT_ASSERT(std::distance(alpha.first, alpha.second) == 1);
  CPPUNIT_ASSERT(std::distance(empty.first, empty.second) == 0);
  CPPUNIT_ASSERT(std::distance(haloed.first, haloed.second) == 1);

  CPPUNIT_ASSERT(helpers::is_zero(omega.first->first));
  CPPUNIT_ASSERT(omega.first->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(
      omega.first->second.nodeIndex == 0 or omega.first->second.nodeIndex == 1
  );
  CPPUNIT_ASSERT(omega.first->second.isNearBorder == 0);
  DivideConquer<CellReference>::const_iterator other = omega.first; ++other;
  CPPUNIT_ASSERT(helpers::is_zero(other->first));
  CPPUNIT_ASSERT(other->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(other->second.nodeIndex == 0 or other->second.nodeIndex == 1);
  CPPUNIT_ASSERT(other->second.nodeIndex != omega.first->second.nodeIndex);
  CPPUNIT_ASSERT(other->second.isNearBorder == 0);

  CPPUNIT_ASSERT(alpha.first->first == LatticeVector(2, 3, -2));
  CPPUNIT_ASSERT(alpha.first->second.nodeIndex == 2);
  CPPUNIT_ASSERT(alpha.first->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(alpha.first->second.isNearBorder == 0);

  CPPUNIT_ASSERT(haloed.first->first == LatticeVector(1, 0, 0));
  CPPUNIT_ASSERT(haloed.first->second.nodeIndex == 3);
  CPPUNIT_ASSERT(haloed.first->second.cellIndex == cellIndex);
  CPPUNIT_ASSERT(haloed.first->second.isNearBorder == CellReference::NORTH);
}

void check_cell(
    DivideConquerCells const &_dnc,
    LatticeVector const &_key,
    site_t _cellIndex,
    site_t _nbNodes
) {
  DivideConquerCells::const_range omega = _dnc(_key);
  CPPUNIT_ASSERT(std::distance(omega.first, omega.second) == _nbNodes);

  std::set<int> nodes;
  for(; omega.first != omega.second; ++omega.first) {
    CPPUNIT_ASSERT(omega.first.GetCellReference().cellIndex == _cellIndex);
    CPPUNIT_ASSERT(
        omega.first.GetCellReference().nodeIndex >= 0
        and omega.first.GetCellReference().nodeIndex <= _nbNodes
    );
    CPPUNIT_ASSERT(
        nodes.count(omega.first.GetCellReference().nodeIndex) == 0);
    nodes.insert(omega.first.GetCellReference().nodeIndex);
  }
  CPPUNIT_ASSERT(nodes.size() == _nbNodes);
}

boost::shared_ptr<CellContainer> CellCellInteractionTests :: fixture(
    PhysicalDistance _cutoff) const {
  boost::shared_ptr<CellContainer> cells(new CellContainer);
  Mesh pancake = pancakeSamosa();
  pancake += LatticePosition(1, 1, 1) * _cutoff * 0.5;
  // safer to clone so cells has its own copy
  cells->push_back(pancake.clone());
  pancake += LatticePosition(3, 0, 1) * _cutoff;
  cells->push_back(pancake.clone());

  return cells;
}

void CellCellInteractionTests :: testAddMeshes() {
  PhysicalDistance const cutoff = 5.0;
  PhysicalDistance const halo = 2.0;
  boost::shared_ptr<CellContainer> cells(fixture(cutoff));
  Mesh pancake = pancakeSamosa();

  DivideConquerCells dnc(cells, cutoff, halo);
  check_cell(dnc, LatticeVector(0, 0, 0), 0, pancake.GetVertices().size());
  check_cell(dnc, LatticeVector(3, 0, 1), 1, pancake.GetVertices().size());
}

void CellCellInteractionTests :: testIterator() {
  PhysicalDistance const cutoff = 5.0;
  PhysicalDistance const halo = 2.0;
  boost::shared_ptr<CellContainer> cells(fixture(cutoff));
  DivideConquerCells dnc(cells, cutoff, halo);

  std::set<LatticePosition const*> allnodes;
  typedef MeshData::t_Vertices::const_iterator vertex_iterator;
  vertex_iterator i_vert = cells->begin()->GetVertices().begin();
  vertex_iterator i_vertend = cells->begin()->GetVertices().end();
  for(; i_vert != i_vertend; ++i_vert)
    allnodes.insert(&(*i_vert));
  i_vert = cells->back().GetVertices().begin();
  i_vertend = cells->back().GetVertices().end();
  for(; i_vert != i_vertend; ++i_vert)
    allnodes.insert(&(*i_vert));


  DivideConquerCells::const_iterator i_first = dnc.begin();
  DivideConquerCells::const_iterator const i_end = dnc.end();
  for(; i_first != i_end; ++i_first) {
    CPPUNIT_ASSERT(not helpers::is_zero(*i_first));
    CPPUNIT_ASSERT(allnodes.count(&(*i_first)) == 1);
    allnodes.erase(&(*i_first));
  }
  CPPUNIT_ASSERT(allnodes.size() == 0);
}

void CellCellInteractionTests :: testUpdate() {
  PhysicalDistance const cutoff = 5.0;
  PhysicalDistance const halo = 2.0;
  boost::shared_ptr<CellContainer> cells(fixture(cutoff));
  DivideConquerCells dnc(cells, cutoff, halo);

  LatticeVector const zero(0, 0, 0);
  LatticeVector const notzero(3, 0, 1);
  CPPUNIT_ASSERT(dnc.size() == 6);
  CPPUNIT_ASSERT(std::distance(dnc(zero).first, dnc(zero).second) == 3);
  CPPUNIT_ASSERT(std::distance(dnc(notzero).first, dnc(notzero).second) == 3);

  // Move to new box
  LatticeVector const newbox(-1, 1, 2);
  LatticePosition const center(0.5 * cutoff);
  LatticePosition const offhalo = LatticePosition(newbox) * cutoff + center;
  cells->front().GetVertices().front() = offhalo;

  dnc.update();
  CPPUNIT_ASSERT(dnc.size() == 6);
  CPPUNIT_ASSERT(std::distance(dnc(notzero).first, dnc(notzero).second) == 3);
  CPPUNIT_ASSERT(std::distance(dnc(zero).first, dnc(zero).second) == 2);
  CPPUNIT_ASSERT(std::distance(dnc(newbox).first, dnc(newbox).second) == 1);
  CPPUNIT_ASSERT(helpers::is_zero(*dnc(newbox).first - offhalo));
  CPPUNIT_ASSERT(not dnc(newbox).first.IsNearBorder());

  // Move near boundary
  LatticePosition const inhalo
    = LatticePosition(0.25, 0, 0.25) * cutoff + offhalo;
  cells->front().GetVertices().front() = inhalo;
  dnc.update();
  CPPUNIT_ASSERT(dnc.size() == 6);
  CPPUNIT_ASSERT(std::distance(dnc(notzero).first, dnc(notzero).second) == 3);
  CPPUNIT_ASSERT(std::distance(dnc(zero).first, dnc(zero).second) == 2);
  CPPUNIT_ASSERT(std::distance(dnc(newbox).first, dnc(newbox).second) == 1);
  CPPUNIT_ASSERT(helpers::is_zero(*dnc(newbox).first - inhalo));
  CPPUNIT_ASSERT(dnc(newbox).first.IsNearBorder());
  CPPUNIT_ASSERT(dnc(newbox).first.IsNearBorder(CellReference::TOP));
  CPPUNIT_ASSERT(not dnc(newbox).first.IsNearBorder(CellReference::BOTTOM));
  CPPUNIT_ASSERT(
      dnc(newbox).first.GetNearBorder()
      == (CellReference::TOP | CellReference::EAST)
  );
}


CPPUNIT_TEST_SUITE_REGISTRATION(CellCellInteractionTests);
}}}

#endif // ONCE
