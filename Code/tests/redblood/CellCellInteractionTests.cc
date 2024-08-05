// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/Cell.h"
#include "redblood/CellCell.cc"

#include "tests/helpers/ApproxVector.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::redblood;

    LatticeDistance constexpr cutoff = 5.0;
    LatticeDistance constexpr halo = 2.0;

    TEST_CASE("testBoxHaloTooBig", "[redblood]") {
      DivideConquer<int> dnc(cutoff);
      LatticeVector const key(0, 1, -2);
      LatticePosition const position = (LatticePosition(key) + LatticePosition(0.6, 0.5, 0.5))
	* cutoff;

      REQUIRE(
	      (static_cast<size_t>(Borders::CENTER)
	       bitor static_cast<size_t>(Borders::NORTH)
	       bitor static_cast<size_t>(Borders::SOUTH)
	       bitor static_cast<size_t>(Borders::EAST)
	       bitor static_cast<size_t>(Borders::WEST)
	       bitor static_cast<size_t>(Borders::TOP)
	       bitor static_cast<size_t>(Borders::BOTTOM))
	      == 
	      figureNearness(dnc, position, cutoff * 1.001)
	      );
      REQUIRE(((size_t) Borders::CENTER bitor (size_t)Borders::TOP) ==
	      figureNearness(dnc, position, cutoff * 0.499));
    }

    TEST_CASE("testBoxHalo", "[redblood]") {
      DivideConquer<int> dnc(5.0);
      LatticeVector const key(0, 1, -2);
      LatticePosition const center = (LatticePosition(key) + LatticePosition(0.5, 0.5, 0.5))
	* dnc.GetBoxSize();
      REQUIRE(figureNearness(dnc, center, 2.0) == (size_t) Borders::CENTER);

      for (size_t d(1); d < (1 << 6); d <<= 1)
        {
          LatticePosition const disp = direction < LatticePosition::value_type > (d) * 0.6;
          size_t const nearness = figureNearness(dnc, center + disp, 2.0);
          REQUIRE(nearness == (d bitor (size_t) Borders::CENTER));
        }

      LatticePosition const mult = direction < LatticePosition::value_type
            > (Borders::TOP) * 0.6 + direction < LatticePosition::value_type
            > (Borders::NORTH) * 0.6 + direction < LatticePosition::value_type
            > (Borders::EAST) * 0.6;
      size_t const actual = figureNearness(dnc, center + mult, 2.0);
      size_t const expected = size_t(Borders::TOP) bitor size_t(Borders::NORTH)
	bitor size_t(Borders::EAST) bitor size_t(Borders::CENTER);
      REQUIRE(actual == expected);
    }

    TEST_CASE("testAddNodes", "[redblood]") {
      DivideConquer < CellReference > dnc(cutoff);

      // Adds nodes, last in halo
      LatticePosition const center = LatticePosition(1, 1, 1) * (0.5 * cutoff);
      LatticeDistance const offhalo = (cutoff * 0.5 - 1.1 * halo) / cutoff;
      LatticeDistance const inhalo = (cutoff * 0.5 - 0.9 * halo) / cutoff;
      MeshData::Vertices vertices;
      vertices.push_back(center);
      vertices.push_back(center + LatticePosition(offhalo, 0, 0) * cutoff);
      vertices.push_back(center + LatticePosition(2, offhalo + 3.0, -2) * cutoff);
      vertices.push_back(center + LatticePosition(1, 0, 0) * cutoff + direction
			 < LatticePosition::value_type > (Borders::NORTH) * inhalo * cutoff);

      CellContainer cells;
      std::shared_ptr<Cell> intel(new Cell(Mesh(MeshData())));
      cells.insert(intel);

      initializeCells(dnc, vertices, cells.begin(), halo);
      REQUIRE(dnc.size() == vertices.size());

      DivideConquer<CellReference>::const_range const omega = dnc.equal_range(LatticeVector(0,
											    0,
											    0));
      DivideConquer<CellReference>::const_range const empty = dnc.equal_range(LatticeVector(10,
											    10,
											    10));
      DivideConquer<CellReference>::const_range const alpha = dnc.equal_range(LatticeVector(2,
											    3,
											    -2));
      DivideConquer<CellReference>::const_range const haloed = dnc.equal_range(LatticeVector(1,
											     0,
											     0));

      REQUIRE(std::distance(omega.first, omega.second) == 2);
      REQUIRE(std::distance(alpha.first, alpha.second) == 1);
      REQUIRE(std::distance(empty.first, empty.second) == 0);
      REQUIRE(std::distance(haloed.first, haloed.second) == 1);

      REQUIRE(LatticeVector::Zero() == omega.first->first);
      REQUIRE(omega.first->second.cellIterator == cells.begin());
      REQUIRE((omega.first->second.nodeIndex == 0 or omega.first->second.nodeIndex == 1));
      REQUIRE(omega.first->second.nearBorder == (size_t) Borders::CENTER);
      DivideConquer<CellReference>::const_iterator other = omega.first;
      ++other;
      REQUIRE(LatticeVector::Zero() == other->first);
      REQUIRE(other->second.cellIterator == cells.begin());
      REQUIRE((other->second.nodeIndex == 0 or other->second.nodeIndex == 1));
      REQUIRE(other->second.nodeIndex != omega.first->second.nodeIndex);
      REQUIRE(other->second.nearBorder == (size_t) Borders::CENTER);

      REQUIRE(alpha.first->first == LatticeVector(2, 3, -2));
      REQUIRE(alpha.first->second.nodeIndex == 2);
      REQUIRE(alpha.first->second.cellIterator == cells.begin());
      REQUIRE(alpha.first->second.nearBorder == (size_t) Borders::CENTER);

      REQUIRE(LatticeVector::value_type(1) == haloed.first->first.x());
      REQUIRE(LatticeVector::value_type(0) == haloed.first->first.y());
      REQUIRE(LatticeVector::value_type(0) == haloed.first->first.z());
      REQUIRE(site_t(3) == haloed.first->second.nodeIndex);
      REQUIRE(cells.begin() == haloed.first->second.cellIterator);
      REQUIRE((size_t(Borders::NORTH) bitor size_t(Borders::CENTER))
	      == haloed.first->second.nearBorder);
    }

    void checkCell(DivideConquerCells const &dnc, LatticeVector const &key,
		   CellContainer::value_type cell, site_t nbNodes)
    {
      DivideConquerCells::const_range omega = dnc(key);
      REQUIRE(std::distance(omega.first, omega.second) == nbNodes);

      std::set<int> nodes;

      for (; omega.first != omega.second; ++omega.first) {
	REQUIRE(*omega.first.GetCellReference().cellIterator == cell);
	REQUIRE(omega.first.GetCellReference().nodeIndex >= 0);
	REQUIRE(omega.first.GetCellReference().nodeIndex <= nbNodes);
	REQUIRE(nodes.count(omega.first.GetCellReference().nodeIndex) == 0);
	nodes.insert(omega.first.GetCellReference().nodeIndex);
      }

      REQUIRE(nodes.size() == static_cast<std::size_t>(nbNodes));
    }

    TEST_CASE("testAddMeshes", "[redblood]")
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        Mesh pancake = pancakeSamosa();

        REQUIRE(cells.size() == 2);
        auto first = *cells.begin();
        auto second = *std::next(cells.begin());
        if (first->GetBarycentre().GetMagnitude() > second->GetBarycentre().GetMagnitude())
          std::swap(first, second);

        DivideConquerCells dnc(cells, cutoff, halo);
        checkCell(dnc, LatticeVector(0, 0, 0), first, pancake.GetVertices().size());
        checkCell(dnc, LatticeVector(3, 0, 1), second, pancake.GetVertices().size());
      }
    
    auto const approx_zero = ApproxVector<double>{0}.Margin(1e-8);

      TEST_CASE("testIterator", "[redblood]")
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);

        std::set<LatticePosition const *> allnodes;
        typedef MeshData::Vertices::const_iterator vertex_iterator;
        vertex_iterator i_vert = (*cells.begin())->GetVertices().begin();
        vertex_iterator i_vertend = (*cells.begin())->GetVertices().end();

        for (; i_vert != i_vertend; ++i_vert)
        {
          allnodes.insert(& (*i_vert));
        }

        i_vert = (*std::next(cells.begin()))->GetVertices().begin();
        i_vertend = (*std::next(cells.begin()))->GetVertices().end();

        for (; i_vert != i_vertend; ++i_vert)
        {
          allnodes.insert(& (*i_vert));
        }

        DivideConquerCells::const_iterator i_first = dnc.begin();
        DivideConquerCells::const_iterator const i_end = dnc.end();

        for (; i_first != i_end; ++i_first)
        {
          REQUIRE(approx_zero != *i_first);
          REQUIRE(allnodes.count(& (*i_first)) == 1);
          allnodes.erase(& (*i_first));
        }

        REQUIRE(allnodes.size() == 0);
      }

      TEST_CASE("testUpdate", "[redblood]")
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);
        // Figures out which cell is which
        auto firstCell = *cells.begin();
        auto secondCell = *std::next(cells.begin());
        if (firstCell->GetBarycentre().GetMagnitude() > secondCell->GetBarycentre().GetMagnitude())
        {
          std::swap(firstCell, secondCell);
        }

        LatticeVector const zero(0, 0, 0);
        LatticeVector const notzero(3, 0, 1);
        REQUIRE(dnc.size() == static_cast<unsigned long>(6));
        REQUIRE(3l == std::distance(dnc(zero).first, dnc(zero).second));
        REQUIRE(3l == std::distance(dnc(notzero).first, dnc(notzero).second));

        // Move to new box
        LatticeVector const newbox(-1, 1, 2);
        LatticePosition const center(0.5 * cutoff);
        LatticePosition const offhalo = LatticePosition(newbox) * cutoff + center;
        firstCell->GetVertices().front() = offhalo;

        dnc.update();
        REQUIRE(dnc.size() == static_cast<unsigned long>(6));
        REQUIRE(3l == std::distance(dnc(notzero).first, dnc(notzero).second));
        REQUIRE(2l == std::distance(dnc(zero).first, dnc(zero).second));
        REQUIRE(1l == std::distance(dnc(newbox).first, dnc(newbox).second));
        REQUIRE(approx_zero == *dnc(newbox).first - offhalo);
        REQUIRE(not dnc(newbox).first.IsNearBorder());

        // Move near boundary
        LatticePosition const inhalo = LatticePosition(0.25, 0, 0.25) * cutoff + offhalo;
        firstCell->GetVertices().front() = inhalo;
        dnc.update();
        REQUIRE(dnc.size() == static_cast<unsigned long>(6));
        REQUIRE(std::distance(dnc(notzero).first, dnc(notzero).second) == 3l);
        REQUIRE(std::distance(dnc(zero).first, dnc(zero).second) == 2l);
        REQUIRE(std::distance(dnc(newbox).first, dnc(newbox).second) == 1l);
        REQUIRE(approx_zero == *dnc(newbox).first - inhalo);
        REQUIRE(dnc(newbox).first.IsNearBorder());
        REQUIRE(dnc(newbox).first.IsNearBorder(Borders::TOP));
        REQUIRE(not dnc(newbox).first.IsNearBorder(Borders::BOTTOM));
        REQUIRE(dnc(newbox).first.GetNearBorder() == (size_t(Borders::TOP) bitor size_t(Borders::EAST)
						      bitor size_t(Borders::CENTER)));
      }

    TEST_CASE("Updating lent cell does not keep extra nodes when number of nodes change in lent cells", "[redblood]") {
      auto mesh = redblood::pancakeSamosa();
      auto cell = std::make_shared<redblood::Cell>(mesh);
      CellContainer const empty;
      DivideConquerCells dnc(empty, cutoff, halo);

      parallel::LentCells lent;
      lent[0].insert(cell);
      dnc.update(parallel::ExchangeCells::ChangedCells(empty, empty, lent));
      REQUIRE(cell->GetVertices().size() == dnc.size());

      cell->GetVertices().pop_back();
      cell->GetVertices().push_back({4, 4, 4});
      dnc.update(parallel::ExchangeCells::ChangedCells(empty, empty, lent));
      REQUIRE(cell->GetVertices().size() == dnc.size());

      cell->GetVertices().resize(1);
      dnc.update(parallel::ExchangeCells::ChangedCells(empty, empty, lent));
      REQUIRE(cell->GetVertices().size() == dnc.size());
    }

      TEST_CASE("testPairIteratorNoPairs", "[redblood]")
      {
        auto cells = TwoPancakeSamosas<>(cutoff);
        DivideConquerCells dnc(cells, cutoff, halo);

        // Test when iterating over nothing
        {
          DivideConquerCells::pair_range range(dnc, dnc.end(), dnc.end(), 0.5);
          REQUIRE(not range.is_valid());
        }

        // Test when no pair are within given range
        {
          DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
          REQUIRE(not range.is_valid());
        }
      }

      TEST_CASE("testPairIteratorSameMesh", "[redblood]")
      {
        auto cells = TwoPancakeSamosas<>(cutoff);

        // Move one node closer  to the other
        LatticePosition const n0 = (*cells.begin())->GetVertices()[0];
        LatticePosition const n1 = (*cells.begin())->GetVertices()[1];
        (*cells.begin())->GetVertices()[1] = (n1 - n0).GetNormalised() * 0.3 + n0;
        DivideConquerCells dnc(cells, cutoff, halo);

        DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
        DivideConquerCells::const_iterator iterator(dnc.begin());
        REQUIRE(not range.is_valid());
      }

    TEST_CASE("testPairIteratorSinglePair", "[redblood]")
    {
      // There is only one pair and they are in the same divide and conquer box
      auto cells = TwoPancakeSamosas<>(cutoff);

      // Move one node closer  to the other
      auto const &firstCell = *cells.begin();
      auto const &secondCell = *next(cells.begin());
      LatticePosition const n0 = firstCell->GetVertices()[0];
      LatticePosition const n1 = secondCell->GetVertices()[1];
      secondCell->GetVertices()[1] = (n1 - n0).GetNormalised() * 0.3 + n0;
      DivideConquerCells dnc(cells, cutoff, halo);

      DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
      REQUIRE(range.is_valid());
      LatticePosition const& firstNode = range->first.GetCell() == firstCell ?
	*range->first :
	*range->second;
      LatticePosition const& secondNode = range->first.GetCell() == firstCell ?
	*range->second :
	*range->first;
      REQUIRE(ApproxV(n0) == firstNode);
      REQUIRE(ApproxV(secondCell->GetVertices()[1]) == secondNode);
      REQUIRE(not ++range);
      REQUIRE(not range.is_valid());
    }

    TEST_CASE("testPairIteratorOnePairPerBox", "[redblood]")
    {
      // There three pairs and they are each in different boxe, but each contained
      // within one box
      auto cells = TwoPancakeSamosas<>(cutoff);

      // Only one pair, and each in a separate box
      LatticePosition const n0(2 * cutoff - 0.1, 4.5 * cutoff, 4.5 * cutoff);
      LatticePosition const n1(2 * cutoff + 0.1, 4.5 * cutoff, 4.5 * cutoff);
      auto const &firstCell = *cells.begin();
      auto const &secondCell = *next(cells.begin());
      firstCell->GetVertices().front() = n0;
      secondCell->GetVertices().front() = n1;

      DivideConquerCells dnc(cells, cutoff, halo);
      // Checks that fixture is what I think it is
      LatticeVector const N0(1, 4, 4);
      LatticeVector const N1(2, 4, 4);
      REQUIRE(std::distance(dnc(N0).first, dnc(N0).second) == 1);
      REQUIRE(std::distance(dnc(N1).first, dnc(N1).second) == 1);

      DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
      REQUIRE(range.is_valid());

      auto const firstNode = range->first.GetCell() == firstCell ?
	*range->first :
	*range->second;
      auto const secondNode = range->first.GetCell() == firstCell ?
	*range->second :
	*range->first;
      REQUIRE(ApproxV(n0) == firstNode);
      REQUIRE(ApproxV(n1) == secondNode);
      REQUIRE(not ++range);
      REQUIRE(not range.is_valid());
    }

    TEST_CASE("testPairIteratorDiagonalBoxes", "[redblood]")
    {
      // There three pairs and they are each in different boxe, but each contained
      // within one box
      auto cells = TwoPancakeSamosas<>(cutoff);

      // Only one pair, and each in diagonally separate box
      LatticePosition const n0(2 * cutoff - 0.1, 4 * cutoff - 0.1, 8 * cutoff - 0.1);
      LatticePosition const n1(2 * cutoff + 0.1, 4 * cutoff + 0.1, 8 * cutoff + 0.1);
      auto const &firstCell = *cells.begin();
      auto const &secondCell = *next(cells.begin());
      firstCell->GetVertices().front() = n0;
      secondCell->GetVertices().front() = n1;

      DivideConquerCells dnc(cells, cutoff, halo);
      // Checks that fixture is what I think it is
      LatticeVector const N0(1, 3, 7);
      LatticeVector const N1(2, 4, 8);
      REQUIRE(std::distance(dnc(N0).first, dnc(N0).second) == 1);
      REQUIRE(std::distance(dnc(N1).first, dnc(N1).second) == 1);
      REQUIRE((size_t(Borders::CENTER) bitor size_t(Borders::TOP)
	       bitor size_t(Borders::NORTH) bitor size_t(Borders::EAST))
	      == dnc(N0).first.GetNearBorder());
      REQUIRE((size_t(Borders::CENTER) bitor size_t(Borders::BOTTOM)
	       bitor size_t(Borders::SOUTH) bitor size_t(Borders::WEST))
	      == dnc(N1).first.GetNearBorder());

      DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
      REQUIRE(range.is_valid());

      auto const firstNode = range->first.GetCell() == firstCell ?
	*range->first :
	*range->second;
      auto const secondNode = range->first.GetCell() == firstCell ?
	*range->second :
	*range->first;
      REQUIRE(ApproxV(n0) == firstNode);
      REQUIRE(ApproxV(n1) == secondNode);
      REQUIRE(not ++range);
      REQUIRE(not range.is_valid());
    }

    TEST_CASE("testPairIteratorBoxHalo", "[redblood]")
    {
      // There three pairs and they are each in different boxe, but each contained
      // within one box
      auto cells = TwoPancakeSamosas<>(cutoff);

      // Only one pair, and each in a separate box
      LatticePosition const n0(2 * cutoff - 0.1, 4.5 * cutoff, 4.5 * cutoff);
      LatticePosition const n1(2 * cutoff + 0.1, 4.5 * cutoff, 4.5 * cutoff);
      auto const &firstCell = *cells.begin();
      auto const &secondCell = *next(cells.begin());
      firstCell->GetVertices().front() = n0;
      secondCell->GetVertices().front() = n1;

      DivideConquerCells dnc(cells, cutoff, halo);
      // Checks that fixture is what I think it is
      LatticeVector const N0(1, 4, 4);
      LatticeVector const N1(2, 4, 4);
      REQUIRE(std::distance(dnc(N0).first, dnc(N0).second) == 1);
      REQUIRE(std::distance(dnc(N1).first, dnc(N1).second) == 1);

      DivideConquerCells::pair_range range(dnc, dnc.begin(), dnc.end(), 0.5);
      REQUIRE(range.is_valid());

      auto const firstNode = range->first.GetCell() == firstCell ?
	*range->first :
	*range->second;
      auto const secondNode = range->first.GetCell() == firstCell ?
	*range->second :
	*range->first;
      REQUIRE(ApproxV(n0) == firstNode);
      REQUIRE(ApproxV(n1) == secondNode);
      REQUIRE(not ++range);
      REQUIRE(not range.is_valid());
    }

    TEST_CASE("testAddCell", "[redblood]")
    {
      auto cell = std::make_shared<Cell>(pancakeSamosa());
      *cell += LatticePosition(1, 1, 1) * cutoff * 0.5;

      CellContainer intelWay;
      auto intelPtr = std::make_shared<Cell>(*cell);
      intelWay.insert(intelPtr);
      DivideConquerCells dnc(intelWay, cutoff, halo);
      const LatticePosition zero(0, 0, 0);
      const LatticePosition one(1, 1, 1);
      REQUIRE(3ul == dnc.size());
      REQUIRE(3l == std::distance(dnc(zero).first, dnc(zero).second));

      *cell += LatticePosition(1, 1, 1) * cutoff * 1.5;
      dnc.insert(cell);
      REQUIRE(6ul == dnc.size());
      REQUIRE(3l == std::distance(dnc(zero).first, dnc(zero).second));
      REQUIRE(3l == std::distance(dnc(one).first, dnc(one).second));

      // verifies calling update does nothing
      dnc.update();
      REQUIRE(6ul == dnc.size());
    }

    TEST_CASE("testRemoveCell", "[redblood]")
    {
      const LatticePosition zero(0, 0, 0);
      auto cells = TwoPancakeSamosas<>(cutoff);
      DivideConquerCells dnc(cells, cutoff, halo);
      auto const zeroBox = dnc(zero);
      REQUIRE(3l == std::distance(zeroBox.first, zeroBox.second));
      auto const cell = zeroBox.first.GetCell();
      dnc.remove(cell);
      REQUIRE(3ul == dnc.size());
      REQUIRE(dnc(zero).first == dnc(zero).second);

      // verifies calling update does nothing
      dnc.update();
      REQUIRE(3ul == dnc.size());
    }

    TEST_CASE("testDetailIteratorBase", "[redblood]")
    {
      // Mutable and constant iterators
      using mut = DivideConquerCells::iterator;
      using con = DivideConquerCells::const_iterator;

      auto cells = TwoPancakeSamosas<>(cutoff);
      DivideConquerCells m_dnc(cells, cutoff, halo);
      DivideConquerCells const c_dnc(cells, cutoff, halo);

      auto m_iter = m_dnc.begin();
      auto c_iter = c_dnc.begin();

      static_assert(std::is_same<decltype(m_iter), redblood::detail::iterator_base<DivideConquer<CellReference>>>::value,
		    "");
      static_assert(std::is_same<decltype(c_iter), redblood::detail::iterator_base<DivideConquer<CellReference> const>>::value,
		    "");

      {
	// Copy construction
	mut tmp1{m_iter};
	con tmp2{c_iter};
	con tmp3{m_iter};
	// Not allowed
	// mut tmp4{c_iter};
	REQUIRE(tmp1 == m_iter);
	REQUIRE(tmp2 == c_iter);
      }

      {
	// Copy assignment
	mut tmp1 = m_iter;
	con tmp2 = c_iter;
	// Not allowed
	//mut tmp4 = c_iter;
	REQUIRE(tmp1 == m_iter);
	REQUIRE(tmp2 == c_iter);
      }
    }

  }
}
