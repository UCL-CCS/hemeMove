// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_TESTRESOURCES_MESHES_H
#define HLBGMYTOOL_OCT_TESTRESOURCES_MESHES_H

#include <memory>
#include <vector>

#include <vtkSmartPointer.h>
class vtkPolyData;
class vtkXMLPolyDataReader;

#include "Index.h"
#include "Iolet.h"

#include "TestResources/helpers.h"

namespace hemelb::gmytool::oct {

struct MeshData {
  std::vector<Vector> points;
  std::vector<Index> triangles;
  std::vector<Vector> normals;
  std::vector<int> labels;
  std::vector<Iolet> iolets;
};

class SimpleMeshFactory {
 public:
  static std::shared_ptr<MeshData> MkTrivial();
  static std::shared_ptr<MeshData> MkSphere();
  static std::shared_ptr<MeshData> MkDuct();

 private:
  using PDR_ptr = vtkSmartPointer<vtkXMLPolyDataReader>;
  using PD_ptr = vtkSmartPointer<vtkPolyData>;

  static PD_ptr ReadSphere();

  static std::shared_ptr<MeshData> sphere_mesh;

  static std::shared_ptr<MeshData> VtkToMesh(PD_ptr pd);
};

}  // namespace hemelb::gmytool::oct
#endif
