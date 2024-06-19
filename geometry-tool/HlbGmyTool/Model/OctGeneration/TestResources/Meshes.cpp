// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "TestResources/Meshes.h"

#include <catch2/catch.hpp>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>

namespace hemelb::gmytool::oct {

std::shared_ptr<MeshData> SimpleMeshFactory::MkTrivial() {
  auto ans = std::make_shared<MeshData>();
  ans->points = {
      {1.2, 1.2, 1.2}, {1.2, 1.2, 2.2}, {1.2, 2.2, 1.2}, {1.2, 2.2, 2.2}};
  ans->triangles = {{0, 1, 2}, {2, 1, 3}};
  ans->normals = {{-1, 0, 0}, {-1, 0, 0}};
  ans->labels = {-1, -1};
  return ans;
}

std::shared_ptr<MeshData> SimpleMeshFactory::MkSphere() {
  if (!sphere_mesh) {
    auto vtk = ReadSphere();
    sphere_mesh = VtkToMesh(vtk);
  }
  return sphere_mesh;
}

std::shared_ptr<MeshData> SimpleMeshFactory::MkDuct() {
  auto ans = std::make_shared<MeshData>();
  ans->points = {
      {0.5, 0.5, 2.5},  {0.5, 3.5, 2.5},  {3.5, 3.5, 2.5},  {3.5, 0.5, 2.5},

      {0.5, 0.5, 14.5}, {0.5, 3.5, 14.5}, {3.5, 3.5, 14.5}, {3.5, 0.5, 14.5}};

  ans->triangles = {// -x
                    {0, 4, 1},
                    {4, 5, 1},
                    // +y
                    {1, 5, 2},
                    {5, 6, 2},
                    // +x
                    {2, 6, 3},
                    {6, 7, 3},
                    // -y
                    {3, 7, 0},
                    {7, 4, 0},
                    // inlet
                    {0, 1, 3},
                    {1, 2, 3},
                    // outlet
                    {5, 4, 7},
                    {6, 5, 7}};

  ans->normals = {{-1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, +1, 0},
                  {+1, 0, 0}, {+1, 0, 0}, {0, -1, 0}, {0, -1, 0},
                  {0, 0, -1}, {0, 0, -1}, {0, 0, +1}, {0, 0, +1}};

  ans->labels = {-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 1, 1};
  ans->iolets = {Iolet{0, true}, Iolet{0, false}};
  return ans;
}

auto SimpleMeshFactory::ReadSphere() -> PD_ptr {
  auto reader = PDR_ptr::New();
  const std::string&& sphere = GetResource("sphere.vtp");
  reader->SetFileName(sphere.c_str());
  reader->Update();
  REQUIRE(reader->GetOutput() != nullptr);
  return reader->GetOutput();
}

std::shared_ptr<MeshData> SimpleMeshFactory::VtkToMesh(PD_ptr pd) {
  auto ans = std::make_shared<MeshData>();

  // Points
  auto nPoints = pd->GetNumberOfPoints();
  ans->points.resize(nPoints);
  for (auto i = 0; i < nPoints; ++i) {
    auto pt = pd->GetPoint(i);
    for (auto d = 0; d < 3; ++d)
      ans->points[i][d] = pt[d];
  }

  auto nTri = pd->GetNumberOfCells();
  ans->triangles.resize(nTri);
  auto raw_tris = pd->GetPolys()->GetData();

  ans->normals.resize(nTri);
  auto normals = pd->GetCellData()->GetNormals();

  ans->labels.resize(nTri);
  auto labels = pd->GetCellData()->GetScalars();

  auto cell_idx = 0;
  for (auto i = 0; i < nTri; ++i) {
    REQUIRE(raw_tris->GetTuple1(cell_idx) == 3);
    cell_idx++;

    auto n = normals->GetTuple(i);

    for (auto d = 0; d < 3; ++d) {
      ans->triangles[i][d] = raw_tris->GetTuple1(cell_idx + d);
      ans->normals[i][d] = n[d];
    }
    ans->labels[i] = labels->GetTuple1(i);

    cell_idx += 3;
  }

  return ans;
}

std::shared_ptr<MeshData> SimpleMeshFactory::sphere_mesh = nullptr;
}  // namespace hemelb::gmytool::oct
