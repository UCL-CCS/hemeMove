#include "GeometryGenerator.h"
#include "GeometryWriter.h"

#include "Neighbours.h"
#include "Site.h"
#include "Block.h"
#include "Domain.h"

#include "Debug.h"

#include "io/formats/geometry.h"

#include "vtkPolyDataAlgorithm.h"
#include "vtkOBBTree.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"

using namespace hemelb::io::formats;

GeometryGenerator::GeometryGenerator() :
	ClippedSurface(NULL), IsFirstSite(true) {
	Neighbours::Init();
	this->Locator = vtkOBBTree::New();
	this->Locator->SetTolerance(1e-9);
	this->hitPoints = vtkPoints::New();
	this->hitCellIds = vtkIdList::New();
}

GeometryGenerator::~GeometryGenerator() {
	this->Locator->Delete();
	this->hitPoints->Delete();
	this->hitCellIds->Delete();
}

void GeometryGenerator::Execute() {
	// Build our locator.
	this->Locator->SetDataSet(this->ClippedSurface);
	this->Locator->BuildLocator();

	/*
	 * Get the scalars associated with the surface polygons. The scalars hold
	 * the index of the Iolet which they represent, -1 meaning they aren't an
	 * Iolet.
	 */
	this->IoletIdArray = vtkIntArray::SafeDownCast(this->ClippedSurface->GetCellData()->GetScalars());
	if (this->IoletIdArray == NULL) {
		// TODO: raise an exception
		std::cout << "Error getting Iolet ID array from clipped surface" << std::endl;
	}

	Domain domain(this->VoxelSizeMetres, this->ClippedSurface->GetBounds());

	GeometryWriter writer(this->OutputGeometryFile, domain.GetBlockSize(),
			domain.GetBlockCounts(), domain.GetVoxelSizeMetres(), domain.GetOriginMetres());

	for (BlockIterator blockIt = domain.begin(); blockIt != domain.end(); ++blockIt) {
		// Open the BlockStarted context of the writer; this will
		// deal with flushing the state to the file (or not, in the
		// case where there are no fluid sites).
		BlockWriter* blockWriterPtr = writer.StartNextBlock();
		Block& block = *blockIt;

		for (SiteIterator siteIt = block.begin(); siteIt != block.end(); ++siteIt) {
			Site& site = **siteIt;
			this->ClassifySite(site);

			if (site.IsFluid) {
				blockWriterPtr->IncrementFluidSitesCount();
				WriteFluidSite(*blockWriterPtr, site);
			} else {
				WriteSolidSite(*blockWriterPtr, site);
			}
		}
		blockWriterPtr->Finish();
		delete blockWriterPtr;
	}
	writer.Close();
}

void GeometryGenerator::WriteSolidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int> (geometry::SOLID);
	// That's all in this case.
}

void GeometryGenerator::WriteFluidSite(BlockWriter& blockWriter, Site& site) {
	blockWriter << static_cast<unsigned int> (geometry::FLUID);

	// Iterate over the displacements of the neighbourhood
	for (unsigned int i = 0; i < Neighbours::n; ++i) {
		unsigned int cutType = site.Links[i].Type;

		if (cutType == geometry::CUT_NONE) {
			blockWriter << static_cast<unsigned int> (geometry::CUT_NONE);
		} else if (cutType == geometry::CUT_WALL ||
				   cutType == geometry::CUT_INLET ||
				   cutType == geometry::CUT_OUTLET) {
				blockWriter << static_cast<unsigned int> (cutType);
				if (cutType == geometry::CUT_INLET ||
			        cutType == geometry::CUT_OUTLET) {
				    blockWriter << static_cast<unsigned int> (site.Links[i].IoletId);
				}
				blockWriter << static_cast<float> (site.Links[i].Distance);
		} else {
			// TODO: throw some exception
			std::cout << "Unknown cut type " <<
				static_cast<unsigned int> (cutType) << " for site " <<
				site.GetIndex() << std::endl;
		}
	}
}

bool GeometryGenerator::IsInsideSurface(const Vector& point) {
	if (this->Locator->InsideOrOutside(&point[0]) < 0) {
		// -1 => inside surface
		return true;
	} else {
		return false;
	}
}

bool GeometryGenerator::GetIsFluid(Site& site) {
	if (!site.IsFluidKnown) {
		site.IsFluid = this->IsInsideSurface(site.Position);
		site.IsFluidKnown = true;
	}
	return site.IsFluid;
}

/*
 * Perform classification of the supplied sites. Note that
 * this will alter the connected sites that have yet to be
 * classified, as we wish to examine each link only once.
 *
 * Each site must have its IsFluid and IsEdge flags set, along
 * with the CutDistances array (at appropriate indices),
 * WallDistance/Normal and BoundaryDistance/Normal.
 */
void GeometryGenerator::ClassifySite(Site& site) {

	if (!this->GetIsFluid(site)) {
		// Nothing to do for solid sites
		return;
	}

	// It is fluid, hence we need the vector of links
	site.CreateLinksVector();

	for (NeighbourIterator neighIt = site.beginall(); neighIt != site.endall(); ++neighIt) {
		Site& neigh = *neighIt;
		unsigned int iNeigh = neighIt.GetNeighbourIndex();
		LinkData& link = site.Links[iNeigh];

		if (this->GetIsFluid(neigh)) {
			// Link from fluid to fluid
			link.Type = geometry::CUT_NONE;
		} else {
			// Link to a solid site
			this->Locator->IntersectWithLine(&site.Position[0],
					&neigh.Position[0], this->hitPoints, this->hitCellIds);
			Vector hitPoint;

			// First hit.
			this->hitPoints->GetPoint(0, &hitPoint[0]);

			// This is set in any solid case
			link.Distance = (hitPoint - site.Position).GetMagnitude();
			// The distance is in voxels but must be output as a fraction of
			// the lattice vector. Scale it.
			link.Distance /= Neighbours::norms[iNeigh];

			// The index of the cell in the vtkPolyData that was hit
			int hitCellId = this->hitCellIds->GetId(0);
			// The value associated with that cell, which identifies what was hit.
			int ioletId = this->IoletIdArray->GetValue(hitCellId);

			if (ioletId < 0) {
				// -1 => we hit a wall
				link.Type = geometry::CUT_WALL;
			} else {
				// We hit an inlet or outlet
				Iolet* iolet = this->Iolets[ioletId];
				if (iolet->IsInlet) {
					link.Type = geometry::CUT_INLET;
				} else {
					link.Type = geometry::CUT_OUTLET;
				}
				// Set the Id
				link.IoletId = iolet->Id;
			}
		}
	}
}

