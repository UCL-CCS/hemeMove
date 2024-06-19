// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include "Index.h"

namespace {
using DV = hemelb::util::Vector3D<double>;
}

namespace py = pybind11;

PYBIND11_MODULE(CommonGeneration, mod) {
  mod.doc() = "Common geometry generation code";

  // Map Vector3D exceptions into Python IndexError
  py::register_exception<hemelb::gmytool::IndexError>(mod, "DVIndexError",
                                                      PyExc_IndexError);

  // Expose the required bits of Vector3D<double>
  py::class_<DV>(mod, "DoubleVector")
      // __init__
      .def(py::init(&DV::Zero))                 // Zero
      .def(py::init<double>())                  // All the same
      .def(py::init<double, double, double>())  // Three components
      // Attribute acces to elements
      .def_readwrite("x", &DV::x)
      .def_readwrite("y", &DV::y)
      .def_readwrite("z", &DV::z)
      // get/set via []
      .def("__getitem__", [](DV const& self, int i) { return self[i]; })
      .def("__setitem__", [](DV& self, int i, double v) { self[i] = v; })
      .def(py::self / double());
}
