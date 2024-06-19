// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_OCT_POLYDATAGENERATOR_H
#define HLBGMYTOOL_OCT_POLYDATAGENERATOR_H

#include <vector>
#include "GenerationError.h"

class vtkPolyData;

namespace hemelb::gmytool::oct {
class Iolet;

// Class that manages creating the whole process.
//
// Inputs:
//
// uint16_t NumberOfLevels - number of levels of our octree
//
// uint16_t TriangleLevel - index of the level to sort triangles onto
//
// std::string OutputGeometryFile - path to write the output
//
// std::vector<Iolet> Iolets - info about the inlets and outlets
//
// vtkPolyData* ClippedSurface - the geometry to voxelised.  Must be
// scaled to lattice coordinates such that it lies within a cube of
// size 2**NumberOfLevels with it's origin at (0,0,0).
//
class PolyDataGenerator {
 public:
  PolyDataGenerator();

  void Execute();

  std::string const& GetOutputGeometryFile();
  void SetOutputGeometryFile(std::string const& val);

  void SetIolets(std::vector<Iolet*> const& iv);

  void SetNumberOfLevels(int n);
  void SetTriangleLevel(int n);

  vtkPolyData* GetClippedSurface(void);
  void SetClippedSurface(vtkPolyData* val);

 private:
  // Members set from outside to initialise
  uint16_t NumberOfLevels;
  uint16_t TriangleLevel;
  std::string OutputGeometryFile;
  std::vector<Iolet> Iolets;
  vtkPolyData* ClippedSurface;
};

}  // namespace hemelb::gmytool::oct
#endif
