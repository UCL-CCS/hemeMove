
#include "Site.h"
#include "Block.h"
#include "Domain.h"
#include "Debug.h"

Domain::Domain(double VoxelSize, double SurfaceBounds[6],
		unsigned int BlockSize) :
	BlockSize(BlockSize), VoxelSize(VoxelSize) {
	double min, max, size, extra, siteZero;
	int nSites, nBlocks, remainder, totalBlocks = 1;

	for (unsigned int i = 0; i < 3; ++i) {
		min = SurfaceBounds[2 * i];
		max = SurfaceBounds[2 * i + 1];
		size = max - min;
		// int() truncates, we add 2 to make sure there's enough
		// room for the sites just outside.
		nSites = int(size / VoxelSize) + 2;

		// The extra space.
		// Minus one, since we want the number of links, not "fence posts"
		extra = (nSites - 1) * VoxelSize - size;
		// We want to balance this equally with the placement of
		// the first site.
		siteZero = min - 0.5 * extra;

		nBlocks = nSites / BlockSize;
		remainder = nSites % BlockSize;
		if (remainder)
			++nBlocks;
		this->Origin[i] = siteZero;
		this->BlockCounts[i] = nBlocks;
		this->SiteCounts[i] = nBlocks * BlockSize;
		totalBlocks *= nBlocks;
	}
	this->blocks.resize(totalBlocks);
	Log() << "Domain size " << this->BlockCounts << std::endl;
}
//
//Domain::~Domain() {
//	delete this->blocks;
//}

Vector Domain::CalcPositionFromIndex(const Index& index) const {
	Vector ans(index);
	ans *= this->VoxelSize;
	ans += this->Origin;
	return ans;
}

Block& Domain::GetBlock(const Index& index) {
	int i = this->TranslateIndex(index);
	Block* bp = this->blocks[i];
	// If the block hasn't been created yet, do so.
	if (!bp) {
		bp = this->blocks[i] = new Block(*this, index, this->BlockSize);
	}
	return *bp;
}

Site& Domain::GetSite(const Index& gIndex) {
	Block& block = this->GetBlock(gIndex / this->BlockSize);
	return block.GetGlobalSite(gIndex);
}

BlockIterator Domain::begin() {
	return BlockIterator(*this);
}

BlockIterator Domain::end() {
	return BlockIterator(*this, Index(this->BlockCounts[0], 0, 0));
}

BlockIterator::BlockIterator() :
	domain(NULL), current(0, 0, 0) {
	this->maxima = this->domain->BlockCounts - 1;
}

BlockIterator::BlockIterator(Domain& dom) :
	domain(&dom), current(0, 0, 0) {
	this->maxima = this->domain->BlockCounts - 1;
}

BlockIterator::BlockIterator(Domain& dom, const Index& start) :
	domain(&dom), current(start) {
	this->maxima = this->domain->BlockCounts - 1;
}

BlockIterator::BlockIterator(const BlockIterator& other) :
	domain(other.domain), current(other.current), maxima(other.maxima) {
}
//	BlockIterator::~BlockIterator();

BlockIterator& BlockIterator::operator=(const BlockIterator& other) {
	if (this == &other) {
		return (*this);
	}
	this->domain = other.domain;
	this->current = other.current;
	this->maxima = other.maxima;

	return (*this);
}

BlockIterator& BlockIterator::operator++() {
	// Note it is an error to increment an iterator past it's end, so we don't
	// need to handle that case.
	int pos;
	// Delete any unnecessary blocks
	for (int i = this->current.x - 1; i < this->current.x + 1; ++i) {
		if (i < 0)
			continue;
		if (i == this->current.x && i != this->maxima.x)
			continue;

		for (int j = this->current.y - 1; j < this->current.y + 1; ++j) {
			if (j < 0)
				continue;
			if (j == this->current.y && j != this->maxima.y)
				continue;

			for (int k = this->current.z - 1; k < this->current.z + 1; ++k) {
				if (k < 0)
					continue;
				if (k == this->current.z && k != this->maxima.z)
					continue;

				// This block can no longer be reached from the current or later
				// blocks, so delete, and set pointer to null
				pos = this->domain->TranslateIndex(i, j, k);
				delete this->domain->blocks[pos];
				this->domain->blocks[pos] = NULL;
			}
		}
	}

	// Update the index vector
	this->current.z += 1;
	if (this->current.z == this->domain->BlockCounts.z) {
		this->current.z = 0;

		this->current.y += 1;
		if (this->current.y == this->domain->BlockCounts.y) {
			this->current.y = 0;

			this->current.x += 1;
		}
	}
	return *this;
}

bool BlockIterator::operator==(const BlockIterator& other) const {
	return (other.domain == this->domain) && (other.current == this->current);
}

bool BlockIterator::operator!=(const BlockIterator& other) const {
	return !(*this == other);
}

BlockIterator::reference BlockIterator::operator*() {
	return this->domain->GetBlock(this->current);
}

BlockIterator::pointer BlockIterator::operator->() {
	return &(*(*this));
}
