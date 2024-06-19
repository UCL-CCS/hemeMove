// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vtkPolyData.h>

#include "PybindVTKTypeCaster.h"

#include "Iolet.h"
#include "PolyDataGenerator.h"

namespace py = pybind11;
using namespace hemelb::gmytool::oct;

PYBIND11_MODULE(OctGeneration, mod) {
  mod.doc() = "Wrapped octree geometry generation code";

  // Now the Iolet
  py::class_<Iolet>(mod, "Iolet")
      .def(py::init<>())
      .def_readwrite("Id", &Iolet::Id)
      .def_readwrite("IsInlet", &Iolet::IsInlet);

  // PolyDataGenerator
  py::class_<PolyDataGenerator>(mod, "PolyDataGenerator")
      .def(py::init<>())
      .def("GetClippedSurface", &PolyDataGenerator::GetClippedSurface)
      .def("SetClippedSurface", &PolyDataGenerator::SetClippedSurface)
      .def("GetOutputGeometryFile", &PolyDataGenerator::GetOutputGeometryFile)
      .def("SetOutputGeometryFile", &PolyDataGenerator::SetOutputGeometryFile)
      .def("SetIolets", &PolyDataGenerator::SetIolets)
      .def("SetNumberOfLevels", &PolyDataGenerator::SetNumberOfLevels)
      .def("SetTriangleLevel", &PolyDataGenerator::SetTriangleLevel)
      .def("Execute", &PolyDataGenerator::Execute);
  ;
}
