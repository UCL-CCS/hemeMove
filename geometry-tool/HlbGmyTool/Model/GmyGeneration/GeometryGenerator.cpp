// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "GeometryGenerator.h"
#include "GeometryWriter.h"

#include "Block.h"
#include "Domain.h"
#include "Neighbours.h"
#include "Site.h"

#include "Debug.h"

#include <cassert>
#include "io/formats/geometry.h"

namespace hemelb::gmytool::gmy {

using namespace io::formats;

GeometryGenerator::GeometryGenerator() {
  Neighbours::Init();
}

GeometryGenerator::~GeometryGenerator() {}

void GeometryGenerator::PreExecute() {}

void GeometryGenerator::Execute(bool skipNonIntersectingBlocks) {
  if (skipNonIntersectingBlocks) {
    throw GenerationErrorMessage(
        "Skip non intersecting blocks functionality currently not available. "
        "See ticket #651");
  }

  this->PreExecute();
  double bounds[6];
  this->ComputeBounds(bounds);
  Domain domain(this->SiteCounts, this->BlockSize);

  GeometryWriter writer(this->OutputGeometryFile, this->BlockSize,
                        domain.GetBlockCounts());

  for (BlockIterator blockIt = domain.begin(); blockIt != domain.end();
       ++blockIt) {
    // Open the BlockStarted context of the writer; this will
    // deal with flushing the state to the file (or not, in the
    // case where there are no fluid sites).
    BlockWriter* blockWriterPtr = writer.StartNextBlock();
    Block& block = *blockIt;

    int side = 0;  // represents whether the block is inside (-1) outside (+1)
                   // or undetermined (0)

    if (skipNonIntersectingBlocks) {
      side = this->BlockInsideOrOutsideSurface(block);
    } else {  // don't use the optimisation -- check every site
      side = 0;
    }

    switch (side) {
      case 1:
        // Block is entirely outside the domain.
        // We don't have to do anything.
        break;
      case 0:
        // Block has some surface within it.
        for (SiteIterator siteIt = block.begin(); siteIt != block.end();
             ++siteIt) {
          Site& site = **siteIt;
          this->ClassifySite(site);
          // here we should check site
          if (site.IsFluid) {
            blockWriterPtr->IncrementFluidSitesCount();
            WriteFluidSite(*blockWriterPtr, site);
          } else {
            WriteSolidSite(*blockWriterPtr, site);
          }
        }
        break;
      case -1:
        // Block is entirely inside the domain
        for (SiteIterator siteIt = block.begin(); siteIt != block.end();
             ++siteIt) {
          Site& site = **siteIt;
          site.IsFluidKnown = true;
          site.IsFluid = true;
          site.CreateLinksVector();
          for (unsigned int link_index = 0; link_index < site.Links.size();
               ++link_index) {
            site.Links[link_index].Type = geometry::CutType::NONE;
          }
          blockWriterPtr->IncrementFluidSitesCount();
          WriteFluidSite(*blockWriterPtr, site);
        }
        break;
      default:
        break;
    }
    blockWriterPtr->Finish();
    blockWriterPtr->Write(writer);
    delete blockWriterPtr;
  }
  writer.Close();
}

void GeometryGenerator::WriteSolidSite(BlockWriter& blockWriter, Site& site) {
  blockWriter << static_cast<unsigned int>(geometry::SiteType::SOLID);
  // That's all in this case.
}

void GeometryGenerator::WriteFluidSite(BlockWriter& blockWriter, Site& site) {
  blockWriter << static_cast<unsigned int>(geometry::SiteType::FLUID);

  // Iterate over the displacements of the neighbourhood
  for (unsigned int i = 0; i < Neighbours::n; ++i) {
    auto&& cutType = site.Links[i].Type;

    if (cutType == geometry::CutType::NONE) {
      blockWriter << static_cast<unsigned int>(geometry::CutType::NONE);
    } else if (cutType == geometry::CutType::WALL ||
               cutType == geometry::CutType::INLET ||
               cutType == geometry::CutType::OUTLET) {
      blockWriter << static_cast<unsigned int>(cutType);
      if (cutType == geometry::CutType::INLET ||
          cutType == geometry::CutType::OUTLET) {
        blockWriter << static_cast<unsigned int>(site.Links[i].IoletId);
      }
      blockWriter << static_cast<float>(site.Links[i].Distance);
    } else {
      // TODO: throw some exception
      std::cout << "Unknown cut type " << static_cast<unsigned int>(cutType)
                << " for site " << site.GetIndex() << std::endl;
    }
  }

  // Indicate whether the current fluid site has a wall normal available and
  // if so write it.
  if (site.WallNormalAvailable) {
    blockWriter << static_cast<unsigned int>(
        geometry::WallNormalAvailability::AVAILABLE);
    blockWriter << static_cast<float>(site.WallNormal[0]);
    blockWriter << static_cast<float>(site.WallNormal[1]);
    blockWriter << static_cast<float>(site.WallNormal[2]);
  } else {
    blockWriter << static_cast<unsigned int>(
        geometry::WallNormalAvailability::NOT_AVAILABLE);
  }
}

void GeometryGenerator::ComputeAveragedNormal(Site& site) const {
  site.WallNormalAvailable = false;

  if (site.IsFluid) {
    site.WallNormal = {0.0, 0.0, 0.0};

    // Compute a weighted sum of the wall normals available and normalise it.
    for (unsigned neighId = 0; neighId < Neighbours::n; ++neighId) {
      LinkData& link = site.Links[neighId];
      if (link.Type == geometry::CutType::WALL) {
        assert(link.DistanceInVoxels != 0);
        double weight = 1 / link.DistanceInVoxels;
        site.WallNormal[0] += weight * link.WallNormalAtWallCut[0];
        site.WallNormal[1] += weight * link.WallNormalAtWallCut[1];
        site.WallNormal[2] += weight * link.WallNormalAtWallCut[2];
        site.WallNormalAvailable = true;
      }
    }

    // Avoid dividing by 0
    if (site.WallNormalAvailable) {
      site.WallNormal.Normalise();
    }
  }
}

}  // namespace hemelb::gmytool::gmy
