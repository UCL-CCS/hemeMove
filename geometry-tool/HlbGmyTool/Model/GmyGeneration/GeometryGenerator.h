// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_GEOMETRYGENERATOR_H
#define HLBGMYTOOL_GMY_GEOMETRYGENERATOR_H

#include <string>
#include <vector>

#include "GenerationError.h"
#include "Iolet.h"

namespace hemelb::gmytool::gmy {

class GeometryWriter;
class Site;
class BlockWriter;
class Block;

// Base class for creating geometry files.
//
// Assumption: the input object is in a coordinate system in lattice
// units with the site at (0, 0, 0) at the origin.
class GeometryGenerator {
 public:
  GeometryGenerator();
  virtual ~GeometryGenerator();
  void Execute(bool skipNonIntersectingBlocks);

  inline std::string GetOutputGeometryFile(void) {
    return this->OutputGeometryFile;
  }
  inline void SetOutputGeometryFile(std::string val) {
    this->OutputGeometryFile = val;
  }

  inline std::vector<Iolet>& GetIolets() { return this->Iolets; }
  inline std::vector<Iolet> const& GetIolets() const { return this->Iolets; }
  inline void SetIolets(std::vector<Iolet> iv) { this->Iolets = iv; }

  inline void SetSiteCounts(unsigned x, unsigned y, unsigned z) {
    this->SiteCounts[0] = x;
    this->SiteCounts[1] = y;
    this->SiteCounts[2] = z;
  }

  inline void SetBlockSize(unsigned n) { this->BlockSize = n; }

  /**
   * This method implements the algorithm used to approximate the wall normal at
   * a given fluid site. This is done based on the normal of the triangles
   * intersected by each lattice link and the distance to those intersections.
   *
   * Current implementation does a weighted sum of the wall normals. The weights
   * are the reciprocal of cut distances along each link.
   *
   * @param site Site object with the data required by the algorithm.
   */
  void ComputeAveragedNormal(Site& site) const;

 protected:
  virtual void ComputeBounds(double[]) const = 0;
  virtual void PreExecute(void);
  virtual void ClassifySite(Site& site) = 0;
  // virtual void CreateCGALPolygon(void);
  void WriteSolidSite(BlockWriter& blockWriter, Site& site);
  void WriteFluidSite(BlockWriter& blockWriter, Site& site);
  // Members set from outside to initialise
  unsigned SiteCounts[3];
  unsigned BlockSize = 8;
  std::string OutputGeometryFile;
  std::vector<Iolet> Iolets;
  virtual int BlockInsideOrOutsideSurface(const Block& block) = 0;
};

}  // namespace hemelb::gmytool::gmy
#endif  // HLBGMYTOOL_GMY_GEOMETRYGENERATOR_H
