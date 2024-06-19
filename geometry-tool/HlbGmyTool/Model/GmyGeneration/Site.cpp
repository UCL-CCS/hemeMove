// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <limits>

#include "Block.h"
#include "Domain.h"
#include "Iolet.h"
#include "Neighbours.h"
#include "Site.h"

namespace hemelb::gmytool::gmy {

// C'tor
Site::Site(Block& block, Index& index)
    : IsFluidKnown(false),
      IsFluid(false),
      Position(block.GetDomain().CalcPositionWorkingFromIndex(index)),
      block(block),
      index(index),
      WallNormalAvailable(false) {}

// C'tor with index constructed in-place
Site::Site(Block& block, unsigned int i, unsigned int j, unsigned int k)
    : IsFluidKnown(false),
      IsFluid(false),
      block(block),
      index(i, j, k),
      WallNormalAvailable(false) {
  this->Position =
      this->block.GetDomain().CalcPositionWorkingFromIndex(this->index);
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

NeighbourIteratorBase::NeighbourIteratorBase(Site& site, unsigned int startpos)
    : site(&site), domain(&site.block.domain), i(startpos) {
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
NeighbourIteratorBase::NeighbourIteratorBase(const NeighbourIteratorBase& other)
    : site(other.site), domain(other.domain), i(other.i) {}

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
bool NeighbourIteratorBase::operator==(
    const NeighbourIteratorBase& other) const {
  return (other.site == this->site) && (other.i == this->i);
}
// Inequality
bool NeighbourIteratorBase::operator!=(
    const NeighbourIteratorBase& other) const {
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
    int siteBlockIjk = this->domain->TranslateIndex(this->site->block.index);
    int neighBlockIjk =
        this->domain->TranslateIndex(neighIndex / this->domain->GetBlockSize());
    if (neighBlockIjk != siteBlockIjk) {
      // block is different
      return (neighBlockIjk > siteBlockIjk);
    } else {
      // sites are in the same block
      int neighLocalIjk = this->site->block.TranslateIndex(neighIndex);
      int siteLocalIjk = this->site->block.TranslateIndex(this->site->index);
      return (neighLocalIjk > siteLocalIjk);
    }

  } else {
    return false;
  }
}

}  // namespace hemelb::gmytool::gmy
