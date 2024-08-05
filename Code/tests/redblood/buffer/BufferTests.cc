// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
//#include <cppunit/extensions/HelperMacros.h>
#include "redblood/buffer/Buffer.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    class OpenBuffer : public buffer::Buffer
    {
    public:
      OpenBuffer(std::shared_ptr<Cylinder> cyl, CellContainer const& cells = CellContainer()) :
	buffer::Buffer(cyl, cells)
      {
      }
      OpenBuffer(Cylinder const & cyl, CellContainer const& cells = CellContainer()) :
	buffer::Buffer(cyl, cells)
      {
      }

      CellContainer::value_type GetJustDropped() const
      {
	return justDropped;
      }
      LatticeDistance GetOffset() const
      {
	return offset;
      }
      LatticeDistance GetLastZ() const
      {
	return lastZ;
      }
      void SetOffset(LatticeDistance const &off)
      {
	offset = off;
      }
      void SetInteraction(LatticeDistance const &off)
      {
	interactionRadius = off;
      }
      double GetInteraction() const
      {
	return interactionRadius;
      }
      auto insert(CellContainer::value_type cell) -> decltype(CellContainer().insert(cell))
      {
	return buffer::Buffer::insert(cell);
      }
      void updateOffset()
      {
	return buffer::Buffer::updateOffset();
      }
      CellContainer::value_type drop()
      {
	return buffer::Buffer::drop();
      }
      CellContainer::value_type nearestCell() const
      {
	return buffer::Buffer::nearestCell();
      }
      CellContainer::value_type furthestCell() const
      {
	return buffer::Buffer::furthestCell();
      }
      bool isDroppablePosition(CellContainer::value_type const &candidate) const
      {
	return buffer::Buffer::isDroppablePosition(candidate);
      }
      size_t size() const
      {
	return virtuals.size();
      }
      void clear()
      {
	justDropped.reset();
	virtuals.clear();
      }
      void fillBuffer(size_t n)
      {
	buffer::Buffer::fillBuffer(n);
      }
    };

    TEST_CASE("BufferTests", "[redblood][buffer]") {
      
      std::vector<CellContainer::value_type> fillers;
      LatticePosition currentFillerPosition;
      LatticePosition increment;
      Cylinder cylinder;
      cylinder.normal = LatticePosition(1, 0, 0);
      cylinder.origin = LatticePosition(2, 2, 2);
      cylinder.radius = 1.5;
      cylinder.length = 2.0;
      
      auto buffer = std::make_shared<OpenBuffer>(cylinder);
      std::vector<std::shared_ptr<Cell> > cells;
      {
	auto cell = std::make_shared<Cell>(icoSphere());
	*cell -= cell->GetBarycentre();
	cells.push_back(std::make_shared<Cell>(*cell));
	buffer->insert(cells.back());
	*cell += cylinder.normal * 2.0;
	cells.push_back(std::make_shared<Cell>(*cell));
	buffer->insert(cells.back());
	*cell += cylinder.normal * 4.0;
	cells.push_back(std::make_shared<Cell>(*cell));
	buffer->insert(cells.back());

	increment = cylinder.normal * 0.5;
	cell = cells[0]->clone();
	auto filler = [cell, &increment, &currentFillerPosition, &fillers]() {
	  currentFillerPosition += increment;
	  CellContainer::value_type result(cell->clone().release());
	  *result += currentFillerPosition;
	  fillers.push_back(result);
	  return result;
	};
	buffer->SetCellDistributionFunction(filler);
	buffer->SetInteraction(0.6);
	buffer->SetOffset(0);
      }

      auto approx = Approx(0.0).margin(1e-8);

      // Drops cell closest to drop point and translates it to LB coords
      SECTION("testDrop") {
	// Gets original positions
	// We know the cells are in cells are sorted
	auto get_pos = [](CellContainer::value_type a)
	  {
	    return a->GetBarycentre();
	  };
	std::vector<LatticePosition> positions(cells.size());
	std::transform(cells.begin(), cells.end(), positions.begin(), get_pos);

	// Change offset to make it interesting
	buffer->SetOffset(10);

	// justDropped should be invalid since nothing was dropped
	REQUIRE(not buffer->GetJustDropped());

	// Now make sure drop statements yield expected result
	for (std::size_t i = 0; i < cells.size(); ++i) {
	  auto cell = buffer->drop();
	  auto expected_pos = positions[i] + cylinder.normal * buffer->GetOffset()
	    + cylinder.origin;
	  REQUIRE(cell == cells[i]);
	  REQUIRE(cell->GetBarycentre() == ApproxV(expected_pos));
	  REQUIRE(buffer->GetJustDropped() == cell);
	}

	// drop should throw if there are no cells.
	REQUIRE_THROWS_AS(buffer->drop(), Exception);
      }

      SECTION("testOffsetUpdateNoCell") {
	// Remove all cells
	for (size_t i = 0; i < cells.size(); ++i) {
	  buffer->drop();
	}
	// Offset should be reset to zero
	buffer->SetOffset(10);
	buffer->updateOffset();
	REQUIRE(buffer->GetOffset() == approx(0));
      }

      SECTION("testOffsetUpdateNoDropped") {
	// set first cell to something not zero, to make the test interesting
	*cells.front() += cylinder.normal * 1.0;
	buffer->updateOffset();
	// updating should set offet such that dropping the cell will place it at the origin
	auto const cell = buffer->drop();
	auto const barycentre = cell->GetBarycentre();
	REQUIRE(ApproxV(cylinder.origin) == barycentre);
      }

      // Moving perpendicular to the normal should have no effect
      SECTION("testOffsetUpdateDroppedMovedPerpendicular") {
	auto cell = buffer->drop();
	buffer->updateOffset();
	auto const startOffset = buffer->GetOffset();

	*cell += LatticePosition(0, 1, 0) * 1.50;
	buffer->updateOffset();
	REQUIRE(approx(startOffset) == buffer->GetOffset());

	*cell += LatticePosition(0, 0, 1) * 15.0;
	buffer->updateOffset();
	REQUIRE(approx(startOffset) == buffer->GetOffset());
      }

      // Moving backwards when cell are interacting should keep offset fixed (buffer is
      // immoveable)
      SECTION("testOffsetUpdateDroppedMovedBackward") {
	// Special case only if dropped and next cell interact
	buffer->SetInteraction(15e0);
	buffer->updateOffset();
	auto cell = buffer->drop();
	auto const startOffset = buffer->GetOffset();

	*cell += cylinder.normal * 15.0;
	buffer->updateOffset();
	auto const expected = buffer->GetOffset();
	REQUIRE(startOffset == approx(expected));
      }

      // When interacting, cannot move by more than what the cells moves forward
      // And cannot move next cell beyond drop point
      SECTION("testOffsetUpdateDroppedMovedForward") {
	// Special case only if dropped and next cell interact
	buffer->SetInteraction(15e0);
	buffer->updateOffset();
	auto cell = buffer->drop();
	auto const startOffset = buffer->GetOffset();

	auto const delta = 0.1;
	*cell -= cylinder.normal * delta;
	buffer->updateOffset();
	REQUIRE(approx(startOffset - delta) == buffer->GetOffset());

	*cell -= cylinder.normal * 2.0;
	buffer->updateOffset();
	auto const nextDrop = buffer->drop();
	REQUIRE(approx(cylinder.origin.x()) == nextDrop->GetBarycentre().x());
      }

      SECTION("testOutsideInteractionRange") {
	// Special case only if dropped and next cell interact
	buffer->SetInteraction(1e0);

	// Arrange cell positions first
	* (cells[0]) += cylinder.normal * 0.0 - cells[0]->GetBarycentre();
	* (cells[1]) += cylinder.normal * 5.0 - cells[1]->GetBarycentre();
	* (cells[2]) += cylinder.normal * 10.0 - cells[2]->GetBarycentre();

	buffer->updateOffset();
	auto cell = buffer->drop();
	auto const startOffset = buffer->GetOffset();

	REQUIRE(cell == cells[0]);

	// Move forward the first cell, but not so little that they would interact
	*cell -= cylinder.normal * 10.0;
	buffer->updateOffset();
	REQUIRE(approx(startOffset - 5e0) == buffer->GetOffset());

	// Now drop another cell. Offset is limited by interactions.
	cell = buffer->drop();
	REQUIRE(cell == cells[1]);
	buffer->updateOffset();
	REQUIRE(approx(startOffset - 5e0 - 4e0) == buffer->GetOffset());
      }

      SECTION("testNearAndFar") {
	auto const near = buffer->nearestCell();
	REQUIRE(near.get() == cells.front().get());
	auto const far = buffer->furthestCell();
	REQUIRE(far.get() == cells.back().get());
      }

      SECTION("testIsDroppable") {
	auto cell = buffer->nearestCell();
	*cell -= cell->GetBarycentre();
	// Center of cylinder should pass
	*cell += cylinder.normal * cylinder.length * 0.5;
	buffer->SetOffset(0);
	REQUIRE(buffer->isDroppablePosition(cell));
	// Offset moves cylinder towards vascular system
	// Now cell is outside cylinder and outside vascular system
	buffer->SetOffset(cylinder.length);
	REQUIRE(not buffer->isDroppablePosition(cell));
	// Move cell back into cylinder
	*cell -= cylinder.normal * cylinder.length;
	REQUIRE(buffer->isDroppablePosition(cell));
	// Move cell down into vascular system
	*cell -= cylinder.normal * cylinder.length;
	REQUIRE(buffer->isDroppablePosition(cell));
      }

      SECTION("testFillBufferThrowIfNoFunc") {
	auto buf = std::make_shared<OpenBuffer>(cylinder);
	REQUIRE_THROWS_AS(buf->fillBuffer(0), Exception);
      }

      SECTION("testFillBuffer") {
	currentFillerPosition = LatticePosition(0, 0, 0);
	buffer->clear();
	REQUIRE(0 == static_cast<int>(buffer->size()));
	buffer->fillBuffer(0);
	// We should have enough to fill cylinder.length + interaction of these critters
	int const ncells = static_cast<int>( (std::ceil(cylinder.length
							+ buffer->GetInteraction()) / 0.5));
	REQUIRE(ncells == static_cast<int>(buffer->size()));
	// We can still add more if wanted. However, once full, only adds as many as requested.
	buffer->fillBuffer(0);
	REQUIRE(ncells == static_cast<int>(buffer->size()));
	buffer->fillBuffer(2);
	REQUIRE(ncells + 2 == static_cast<int>(buffer->size()));
      }

      SECTION("testNumberOfRequests") {
      	REQUIRE(site_t(0) == buffer->NumberOfRequests());
	buffer->requestNewCells(1);
	REQUIRE(site_t(1) == buffer->NumberOfRequests());
	buffer->requestNewCells(10);
	REQUIRE(site_t(11) == buffer->NumberOfRequests());
      }

      SECTION("testDropsCellsWhenPossible") {
	currentFillerPosition = LatticePosition(0, 0, 0);
	std::vector<CellContainer::value_type> transfered;
	auto addcell = [&transfered](CellContainer::value_type a)
	  {
	    transfered.push_back(a);
	  };
	int const ncells = static_cast<int>( (std::ceil(cylinder.length
							+ buffer->GetInteraction()) / 0.5));

	// There should be no drop if there are no requests
	buffer->clear();
	REQUIRE(static_cast<int>(buffer->size()) == 0);
	(*buffer)(addcell);
	REQUIRE(size_t(0) == transfered.size());

	// Now are going to drop n cells
	buffer->requestNewCells(ncells - 1);
	(*buffer)(addcell);
	// We expect n cells along given direction were transfered
	REQUIRE(size_t(ncells - 1) == transfered.size());
	for (size_t i(0); i < transfered.size(); ++i)
	  {
	    // those cells should be those created when filling the buffer
	    REQUIRE(transfered[i].get() == fillers[i].get());
	  }
	// we expect buffer should still be full
	REQUIRE(static_cast<int>(buffer->size()) >= 1);
      }

      SECTION("testNoDropsCellsWhenNotPossible") {
	currentFillerPosition = LatticePosition(0, 0, 0);
	std::vector<CellContainer::value_type> transfered;
	auto addcell = [&transfered](CellContainer::value_type a)
	  {
	    transfered.push_back(a);
	  };
	int const ncells = static_cast<int>( (std::ceil(cylinder.length
							+ buffer->GetInteraction()) / 0.5));

	buffer->clear();
	buffer->requestNewCells(ncells);
	(*buffer)(addcell);
	// only dropped n - 1 because outside buffer
	REQUIRE(size_t(ncells - 1) == transfered.size());
	REQUIRE(site_t(1) == buffer->NumberOfRequests());
	REQUIRE(static_cast<int>(buffer->size()) >= 1);
	// now move last cell to make space for another, but not enough
	* (transfered.back()) -= cylinder.normal * 0.500 * buffer->GetInteraction();
	(*buffer)(addcell);
	REQUIRE(size_t(ncells - 1) == transfered.size());
	REQUIRE(site_t(1) == buffer->NumberOfRequests());
	// now move last cell to make space for another
	* (transfered.back()) -= cylinder.normal * 0.501 * buffer->GetInteraction();
	buffer->requestNewCells(1);
	(*buffer)(addcell);
	REQUIRE(size_t(ncells) == transfered.size());
	REQUIRE(site_t(1) == buffer->NumberOfRequests());
	REQUIRE(static_cast<int>(buffer->size()) >= 1);
      }

    }
  }
}

