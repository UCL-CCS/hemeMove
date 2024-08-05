# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# import numpy as np
import os.path
import string
from xml.etree.ElementTree import Element, SubElement, ElementTree

# from .Profile import Profile, metre
from .Iolets import Inlet, Outlet
from .Vector import Vector


class HexFloatFormatter(string.Formatter):
    """Write floats as hexadecimal unless other requested"""

    def format_field(self, value, format_spec):
        if isinstance(value, float):
            if format_spec == "" or format_spec == "x":
                return value.hex()
        elif isinstance(value, Vector):
            return self.format(Vector._FORMAT_STRING, value)

        return format(value, format_spec)


class XmlWriter(object):
    VERSION = 6

    FMT = HexFloatFormatter()

    @classmethod
    def _totext(cls, value):
        return cls.FMT.format("{}", value)

    def __init__(self, profile):
        self.profile = profile
        return

    @staticmethod
    def indent(elem, level=0):
        i = "\n" + level * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for elem in elem:
                XmlWriter.indent(elem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = i
        return

    def Write(self):
        root = Element("hemelbsettings", version=str(self.VERSION))
        self.DoSimulation(root)
        self.DoGeometry(root)
        self.DoIolets(root)
        self.DoInitialConditions(root)

        self.indent(root)
        ElementTree(root).write(
            self.profile.OutputXmlFile, encoding="US-ASCII", xml_declaration=True
        )

    def DoSimulation(self, root):
        sim = SubElement(root, "simulation")
        self.QuantityElement(sim, "step_length", self.profile.TimeStepSeconds, "s")
        self.QuantityElement(
            sim,
            "steps",
            round(self.profile.DurationSeconds / self.profile.TimeStepSeconds),
            "lattice",
        )
        self.QuantityElement(sim, "voxel_size", self.profile.VoxelSizeMetres, "m")
        self.QuantityElement(sim, "origin", self.profile.OriginMetres, "m")
        return

    def DoGeometry(self, root):
        geom = SubElement(root, "geometry")
        data = SubElement(geom, "datafile")
        data.set(
            "path",
            os.path.relpath(
                self.profile.OutputGeometryFile,
                os.path.split(self.profile.OutputXmlFile)[0],
            ),
        )
        return geom

    def DoInitialConditions(self, root):
        ic = SubElement(root, "initialconditions")
        pressure = SubElement(ic, "pressure")
        self.QuantityElement(pressure, "uniform", 0.0, "Pa")
        return

    def DoProperties(self, root):
        SubElement(root, "properties")
        return

    def DoIolets(self, root):
        inlets = SubElement(root, "inlets")
        outlets = SubElement(root, "outlets")

        for io in self.profile.Iolets:
            if isinstance(io, Inlet):
                iolet = SubElement(inlets, "inlet")
            elif isinstance(io, Outlet):
                iolet = SubElement(outlets, "outlet")
            else:
                continue

            # At the moment the GUI only allows specification of cosine pressure
            # iolets. Add code here if other types (pressure file, velocity
            # poiseuille/womersley/file) become supported.
            condition = SubElement(
                iolet, "condition", type="pressure", subtype="cosine"
            )
            self.QuantityElement(condition, "mean", io.Pressure.x, "Pa")
            self.QuantityElement(condition, "amplitude", io.Pressure.y, "Pa")
            self.QuantityElement(condition, "phase", io.Pressure.z, "rad")
            self.QuantityElement(condition, "period", 1, "s")

            self.QuantityElement(iolet, "normal", io.Normal, "dimensionless")

            # Scale the centre to metres
            centre = Vector()
            for dim in ("x", "y", "z"):
                setattr(
                    centre,
                    dim,
                    getattr(io.Centre, dim) * self.profile.StlFileUnit.SizeInMetres,
                )

            self.QuantityElement(iolet, "position", centre, "m")
            continue
        return

    def QuantityElement(self, parent, name, value, units):
        return SubElement(parent, name, value=self._totext(value), units=units)
