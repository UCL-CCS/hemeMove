// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "PolyDataGenerator.h"

#include <chrono>
#include "range.hpp"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>

#include "FloodFill.h"
#include "Index.h"
#include "SectionTree.h"
#include "SectionTreeBuilder.h"
#include "SurfaceVoxeliser.h"
#include "TriangleSorter.h"

namespace hemelb::gmytool::oct {

PolyDataGenerator::PolyDataGenerator()
    : NumberOfLevels(0), ClippedSurface(nullptr) {}

std::string const& PolyDataGenerator::GetOutputGeometryFile() {
  return this->OutputGeometryFile;
}
void PolyDataGenerator::SetOutputGeometryFile(std::string const& val) {
  this->OutputGeometryFile = val;
}

void PolyDataGenerator::SetIolets(std::vector<Iolet*> const& iv) {
  this->Iolets.clear();
  for (auto&& ioptr : iv) {
    this->Iolets.push_back(*ioptr);
  }
}

void PolyDataGenerator::SetNumberOfLevels(int n) {
  this->NumberOfLevels = n;
}
void PolyDataGenerator::SetTriangleLevel(int n) {
  this->TriangleLevel = n;
}

vtkPolyData* PolyDataGenerator::GetClippedSurface(void) {
  return this->ClippedSurface;
}
void PolyDataGenerator::SetClippedSurface(vtkPolyData* val) {
  this->ClippedSurface = val;
}

class Timer {
 public:
  typedef std::chrono::high_resolution_clock clock;
  Timer(const char* msg) : message(msg), t0(clock::now()) {}

  float GetSeconds() {
    std::chrono::duration<float> delta_t(clock::now() - t0);
    return delta_t.count();
  }

  template <class OutStream>
  void Report(OutStream& os) {
    os << message << ": " << GetSeconds() << " s" << std::endl;
  }

 private:
  std::string message;
  clock::time_point t0;
};

void PolyDataGenerator::Execute() {
  std::cout << "Begin Execute (times in seconds)" << std::endl;

  Timer total_t("Total");
  Timer convert_t("Convert VTK input data");

  auto npts = this->ClippedSurface->GetNumberOfPoints();
  std::vector<Vector> points;
  points.reserve(npts);
  for (auto i : range(npts)) {
    const auto pt = this->ClippedSurface->GetPoint(i);
    points.emplace_back(pt[0], pt[1], pt[2]);
  }

  auto nTri = this->ClippedSurface->GetNumberOfCells();

  std::vector<Index> triangles;
  std::vector<Vector> normals;
  std::vector<int> labels;

  triangles.reserve(nTri);
  auto raw_tris = this->ClippedSurface->GetPolys()->GetData();
  normals.reserve(nTri);
  auto raw_normals = this->ClippedSurface->GetCellData()->GetNormals();
  labels.reserve(nTri);
  auto raw_labels = this->ClippedSurface->GetCellData()->GetScalars();

  auto cell_idx = 0;
  for (auto i = 0; i < nTri; ++i) {
    cell_idx++;

    triangles.emplace_back(raw_tris->GetTuple1(cell_idx + 0),
                           raw_tris->GetTuple1(cell_idx + 1),
                           raw_tris->GetTuple1(cell_idx + 2));

    auto n = raw_normals->GetTuple(i);
    normals.emplace_back(n[0], n[1], n[2]);
    labels.emplace_back(raw_labels->GetTuple1(i));

    cell_idx += 3;
  }

  convert_t.Report(std::cout);

  Timer tri_sort_t("Sort triangles onto tree");
  auto tree = TrianglesToTreeParallel(this->NumberOfLevels, this->TriangleLevel,
                                      points, triangles, 0);
  tri_sort_t.Report(std::cout);

  Timer voxing_t("Voxelise the surface");
  SurfaceVoxeliser voxer(1 << this->TriangleLevel, points, triangles, normals,
                         labels, this->Iolets);
  auto fluid_tree = voxer(tree, this->TriangleLevel);
  voxing_t.Report(std::cout);

  Timer fill_t("Flood fill");
  // Fill the thing
  FloodFill ff(fluid_tree);
  auto mask_tree = ff();
  fill_t.Report(std::cout);

  Timer section_t("Build the section tree");
  SectionTreeBuilder builder(mask_tree, fluid_tree);
  auto section_tree = builder();
  section_t.Report(std::cout);

  Timer write_t("Write the data");
  section_tree->Write(this->OutputGeometryFile);
  write_t.Report(std::cout);

  total_t.Report(std::cout);
}

}  // namespace hemelb::gmytool::oct
