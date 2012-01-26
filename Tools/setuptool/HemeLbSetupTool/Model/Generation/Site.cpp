#include <limits>

#include "Site.h"
#include "Block.h"
#include "Domain.h"
#include "Iolet.h"
#include "Neighbours.h"

#include "geometry/SiteData.h"

// C'tor
Site::Site(Block& block, Index& index) :
	IsFluidKnown(false), IsFluid(false), IsEdge(false),
			Position(block.GetDomain().CalcPositionFromIndex(index)),
			block(block), index(index) {
	this->Init();
}

// C'tor with index constructed in-place
Site::Site(Block& block, unsigned int i, unsigned int j, unsigned int k) :
	IsFluidKnown(false), IsFluid(false), IsEdge(false), block(block),
			index(i, j, k) {
	this->Position = this->block.GetDomain().CalcPositionFromIndex(this->index);
	this->Init();
}

// private helper for common parts of c'tors
void Site::Init() {
	// Make sure this is false
	this->AdjacentIolet = NULL;
	// We need to find the closest.
	this->BoundaryDistance = std::numeric_limits<double>::infinity();
	this->WallDistance = std::numeric_limits<double>::infinity();

	this->CutDistances.resize(Neighbours::n,
			std::numeric_limits<double>::infinity());
	this->CutCellIds.resize(Neighbours::n, -1);
}

// Compute the type of the site from it's properties
unsigned int Site::GetType() const {
	if (this->IsFluid) {
		if (this->AdjacentIolet) {
			if (this->AdjacentIolet->IsInlet)
				return hemelb::geometry::INLET_TYPE;
			else
				return hemelb::geometry::OUTLET_TYPE;
		} else {
			return hemelb::geometry::FLUID_TYPE;
		}
	} else {
		return hemelb::geometry::SOLID_TYPE;
	}
}

// Compute the full config unsigned int
uint64_t Site::GetConfig() {
	unsigned int type = this->GetType();
	uint64_t cfg = type;
	// If solid, we're done
	if (!this->IsFluid)
		return cfg;

	// Fluid sites now
	if (type == hemelb::geometry::FLUID_TYPE && !this->IsEdge)
		// Simple fluid sites
		return cfg;

	// A complex one

	if (this->IsEdge) {
		// Bit fiddle the boundary config. See comment below
		// by BOUNDARY_CONFIG_MASK for definition.
		unsigned int boundary = 0;
		for (NeighbourIterator neighIt = this->beginall(); neighIt
				!= this->endall(); ++neighIt) {
			if (!neighIt->IsFluid) {
				// If the lattice vector is cut, set the flag
				boundary |= 1 << neighIt.GetNeighbourIndex();
			}
		}
		// Shift the boundary bit field to the appropriate
		// place and set these bits in the type
		cfg |= (uint64_t) boundary << hemelb::geometry::SiteData::BOUNDARY_CONFIG_SHIFT;

		if (boundary)
			// Set this bit if we've hit any solid sites
			cfg |= hemelb::geometry::SiteData::PRESSURE_EDGE_MASK;

	}

	if (type != hemelb::geometry::FLUID_TYPE) {
		// It must be an inlet or outlet
		// Shift the index left and set the bits
		cfg |= (uint64_t) this->BoundaryId << hemelb::geometry::SiteData::BOUNDARY_ID_SHIFT;
	}

	return cfg;
}

const Index Site::GetDomainBlockCount() {
	return block.GetDomain().GetBlockCounts();
}

const int Site::GetDomainBlockSize() {
	return block.GetDomain().GetBlockSize();
}

// Get start and end LaterNeighbourIterators for this
LaterNeighbourIterator Site::begin() {
	return LaterNeighbourIterator(*this);
}
LaterNeighbourIterator Site::end() {
	return LaterNeighbourIterator(*this, Neighbours::n);
}
// Get start and end NeighbourIterators for this
NeighbourIterator Site::beginall() {
	return NeighbourIterator(*this);
}

NeighbourIterator Site::endall() {
	return NeighbourIterator(*this, Neighbours::n);
}

/*
 * Implement the NeighbourIterators
 */

NeighbourIteratorBase::NeighbourIteratorBase(Site& site, unsigned int startpos) :
	site(&site), domain(&site.block.domain), i(startpos) {
	/* Should advance to the first valid neighbour
	 * HOWEVER, we can't do that here, in the abstract class's c'tor
	 * since the subclass c'tors haven't yet executed. The virtual
	 * method table may not be in a well defined state and we need that for
	 * IsCurrent Valid.
	 */
}
void NeighbourIteratorBase::AdvanceToValid() {
	// Advance zero or more times until we're at a vector that is valid.
	while (this->i < Neighbours::n && !IsCurrentValid()) {
		++this->i;
	}
}

// copy c'tor
NeighbourIteratorBase::NeighbourIteratorBase(const NeighbourIteratorBase& other) :
	site(other.site), domain(other.domain), i(other.i) {
}

// assignment operator
NeighbourIteratorBase& NeighbourIteratorBase::operator=(
		const NeighbourIteratorBase& other) {
	if (this == &other) {
		return (*this);
	}
	this->site = other.site;
	this->domain = other.domain;
	this->i = other.i;

	return (*this);
}

// Is the current a site within the full domain?
bool NeighbourIteratorBase::IsCurrentInDomain() {
	Index index = this->site->index + this->GetVector();

	for (unsigned int j = 0; j < 3; ++j) {
		// If ind is out of the Domain, the current is invalid
		if (index[j] < 0 or index[j] >= this->domain->SiteCounts[j])
			return false;
	}
	return true;
}

NeighbourIteratorBase& NeighbourIteratorBase::operator++() {
	// Note it is an error to increment an iterator past it's end, so we don't
	// need to handle that case.

	// Go to the next, then keep advancing if it's not valid
	++this->i;
	this->AdvanceToValid();
	return *this;
}

// Test for equality with another
bool NeighbourIteratorBase::operator==(const NeighbourIteratorBase& other) const {
	return (other.site == this->site) && (other.i == this->i);
}
// Inequality
bool NeighbourIteratorBase::operator!=(const NeighbourIteratorBase& other) const {
	return !(*this == other);
}

// Dereference
NeighbourIteratorBase::reference NeighbourIteratorBase::operator*() {
	return this->domain->GetSite(this->site->index + GetVector());
}
// Member lookup
NeighbourIteratorBase::pointer NeighbourIteratorBase::operator->() {
	return &(*(*this));
}

// Index of neighbour
unsigned int NeighbourIteratorBase::GetNeighbourIndex() {
	return this->i;
}

// Get the lattice vector for the current neighbour
Index NeighbourIteratorBase::GetVector() {
	return Neighbours::vectors[this->i];
}

// Normal iterator just checks the neighbour is in the domain
bool NeighbourIterator::IsCurrentValid() {
	return this->IsCurrentInDomain();
}

// Later one checks
bool LaterNeighbourIterator::IsCurrentValid() {
	Index neighIndex = this->site->index + GetVector();

	if (this->IsCurrentInDomain()) {
		// neighbour's in the domain, check it's block
		int siteBlockIjk =
				this->domain->TranslateIndex(this->site->block.index);
		int neighBlockIjk = this->domain->TranslateIndex(
				neighIndex / this->domain->GetBlockSize());
		if (neighBlockIjk != siteBlockIjk) {
			// block is different
			return (neighBlockIjk > siteBlockIjk);
		} else {
			// sites are in the same block
			int neighLocalIjk = this->site->block.TranslateIndex(neighIndex);
			int siteLocalIjk = this->site->block.TranslateIndex(
					this->site->index);
			return (neighLocalIjk > siteLocalIjk);
		}

	} else {
		return false;
	}
}
