// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vtkPolyData.h>

#include "PybindVTKTypeCaster.h"

#include "CylinderGenerator.h"
#include "GeometryGenerator.h"
#include "Iolet.h"
#include "PolyDataGenerator.h"

namespace py = pybind11;
using namespace hemelb::gmytool::gmy;

PYBIND11_MODULE(GmyGeneration, mod) {
  mod.doc() = "Wrapped GMY geometry generation code";

  // Now the Iolet
  py::class_<Iolet>(mod, "Iolet")
      .def(py::init<>())
      .def_readwrite("Centre", &Iolet::Centre)
      .def_readwrite("Normal", &Iolet::Normal)
      .def_readwrite("Radius", &Iolet::Radius)
      .def_readwrite("Id", &Iolet::Id)
      .def_readwrite("IsInlet", &Iolet::IsInlet);

  // ABC GeometryGenerator
  auto gmy_gen =
      py::class_<GeometryGenerator>(mod, "GeometryGenerator")
          .def("GetOutputGeometryFile",
               &GeometryGenerator::GetOutputGeometryFile)
          .def("SetOutputGeometryFile",
               &GeometryGenerator::SetOutputGeometryFile)
          .def("GetIolets", py::overload_cast<>(&GeometryGenerator::GetIolets))
          .def("SetIolets", &GeometryGenerator::SetIolets)
          .def("SetSiteCounts", &GeometryGenerator::SetSiteCounts)
          .def("SetBlockSize", &GeometryGenerator::SetBlockSize)
          .def("Execute", &GeometryGenerator::Execute);

  // Concrete PolyDataGenerator
  py::class_<PolyDataGenerator>(mod, "PolyDataGenerator", gmy_gen)
      .def(py::init<>())
      .def("GetClippedSurface", &PolyDataGenerator::GetClippedSurface)
      .def("SetClippedSurface", &PolyDataGenerator::SetClippedSurface);

  // Cylinder
  py::class_<CylinderGenerator>(mod, "CylinderGenerator", gmy_gen)
      .def(py::init<>())
      .def("SetCylinderCentre", &CylinderGenerator::SetCylinderCentre)
      .def("SetCylinderAxis", &CylinderGenerator::SetCylinderAxis)
      .def("SetCylinderRadius", &CylinderGenerator::SetCylinderRadius)
      .def("SetCylinderLength", &CylinderGenerator::SetCylinderLength);
}
