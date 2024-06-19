# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from abc import ABC, abstractmethod
import math
import os.path
import sys
import time
from warnings import warn

import numpy as np

np.seterr(divide="ignore")

from vtk import (
    vtkClipPolyData,
    vtkAppendPolyData,
    vtkPlane,
    vtkStripper,
    vtkFeatureEdges,
    vtkPolyDataConnectivityFilter,
    vtkProgrammableFilter,
    vtkTriangleFilter,
    vtkCleanPolyData,
    vtkIntArray,
    vtkPoints,
    vtkPolyData,
    vtkCellArray,
    vtkTransform,
    vtkTransformFilter,
    vtkIdList,
    vtkPolyLine,
    vtkXMLPolyDataWriter,
    vtkAlgorithm,
    vtkImplicitBoolean,
    vtkSphere,
    vtkPolyDataNormals,
    vtkSTLWriter,
)

from vmtk.vtkvmtk import vtkvmtkPolyDataBoundaryExtractor
from vmtk.vtkvmtk import vtkvmtkBoundaryReferenceSystems

from .Iolets import Inlet, Outlet, Iolet
from .Vector import Vector
from .Profile import Profile, metre
from .XmlWriter import XmlWriter

from . import CommonGeneration

# Add Pythonic printing
class DoubleVector(CommonGeneration.DoubleVector):
    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"

    def __repr__(self):
        return f"DoubleVector({self.x}, {self.y}, {self.z})"


def DVfromV(v):
    """Translate a Model.Vector.Vector to a Generation.DoubleVector."""
    return DoubleVector(v.x, v.y, v.z)


try:
    from . import GmyGeneration
except ImportError:
    warn("GMY generation unavailable")
    GmyGeneration = None

try:
    from . import OctGeneration
except ImportError:
    warn("OCT generation unavailable")
    OctGeneration = None

if OctGeneration is None and GmyGeneration is None:
    raise ImportError("No generation extensions available!")


# We have several concrete generators we want, shown in this table:
#
#          | GMY | OCT
# ---------+-----+-----
# PolyData |  X  |  X
# Duct     |  X  |  -
# Cylinder |  X  |  -
#
# (X == implemented, - == not implemented)
#
# The compiled modules for one or other of the back ends might be unavailable
#


# We have four coordinate systems to deal with.
#
# 1. The physical system in metres.
#
# 2. The profile / surface system, which is a scaling of the physical
#    system (by profile.StlFileUnit.SizeInMetres)
#
# 3. The scaled system which is a simple scaling to lattice units.
#
# 4. The working system which is (2) shifted such that a site's lattice
#    indices (i,j,k) are the same as it's coordinates


class GeometryGenerator(ABC):
    def __getattr__(self, attr):
        """Delegate unknown attribute access to profile object."""
        return getattr(self._profile, attr)

    def _MakeIoletProxies(self):
        # Construct the Iolet structs
        nIn = 0
        nOut = 0
        ioletProxies = []
        for io in self.Iolets:
            proxy = self.genext.Iolet()

            self._AddIoletProperties(io, proxy)

            if isinstance(io, Inlet):
                io.Id = proxy.Id = nIn
                proxy.IsInlet = True
                nIn += 1
            elif isinstance(io, Outlet):
                io.Id = proxy.Id = nOut
                proxy.IsInlet = False
                nOut += 1

            ioletProxies.append(proxy)
        return ioletProxies

    @abstractmethod
    def _AddIoletProperties(self, io, proxy):
        pass

    def _SetGeneratorProperties(self):
        """Called at the start of execute.

        Sets up the C++ generator based on our state.
        """
        self.generator.SetOutputGeometryFile(str(self.GeneratorOutputFile))
        self.generator.SetIolets(self._MakeIoletProxies())
        self._SetExtensionCommonProperties()
        self._SetGeometryTypeProperties()

    @abstractmethod
    def _SetExtensionCommonProperties(self):
        pass

    @abstractmethod
    def _SetGeometryTypeProperties(self):
        pass

    @property
    @abstractmethod
    def GeneratorOutputFile(self):
        """Because the file that the C++ writes might not be .gmy"""
        pass

    @abstractmethod
    def _Preprocess(self):
        """Do any preprocessing required (e.g. clip polydata)."""
        pass

    def Execute(self):
        """Pass our state to C++ implementation and run."""
        self._Preprocess()
        self._SetGeneratorProperties()
        t = Timer()
        t.Start()
        self._DoExecute()
        t.Stop()
        print("Setup time: %f s" % t.GetTime())
        mem, unit = t.GetMem()
        if mem is not None:
            print(f"Memory used by generator: {mem} {unit}")
        XmlWriter(self).Write()

    @abstractmethod
    def _DoExecute(self):
        """Different call signatures in C++"""
        pass


class GmyMixin:
    """Behaviour for the direct GMY generation extension."""

    def __init__(self, *args, **kwargs):
        if GmyGeneration is None:
            raise RuntimeError("The geometry generator module is not available")
        self.genext = GmyGeneration

        self.skipNonIntersectingBlocks = False
        self.BlockSize = 8
        super().__init__(*args, **kwargs)

    @property
    def GeneratorOutputFile(self):
        return self._profile.OutputGeometryFile

    def _SetExtensionCommonProperties(self):
        self.generator.SetSiteCounts(*self.SiteCounts)
        self.generator.SetBlockSize(self.BlockSize)

    def _AddIoletProperties(self, io, proxy):
        proxy.Centre = DVfromV(io.Centre) / self.VoxelSize
        proxy.Normal = DVfromV(io.Normal)
        proxy.Radius = io.Radius / self.VoxelSize

    def _DoExecute(self):
        self.generator.Execute(self.skipNonIntersectingBlocks)

    def _RoundUpSites(self, nsites):
        return ((nsites - 1) // self.BlockSize + 1) * self.BlockSize


class OctMixin:
    """Behaviour for the octree generation extension."""

    def __init__(self, *args, **kwargs):
        if OctGeneration is None:
            raise RuntimeError("The Octree generator module is not available")
        self.genext = OctGeneration
        super().__init__(*args, **kwargs)

    @property
    def OutputOctFile(self):
        base, ext = os.path.splitext(self._profile.OutputGeometryFile)
        return base + ".oct"

    @property
    def GeneratorOutputFile(self):
        return self.OutputOctFile

    def _SetExtensionCommonProperties(self):
        pass

    def _DoExecute(self):
        self.generator.Execute()

    def _RoundUpSites(self, nsites):
        # Get the biggest of the sides
        max_sites = int(nSites.max())
        # Now round this up to the next power of two
        nBits = max_sites.bit_length()
        if max_sites > 2 ** (nBits - 1):
            max_sites = 2**nBits
            pass
        self.NumberOfLevels = nBits

        self.CubeSize = max_sites


class PolyDataGenerator(GeometryGenerator):
    def __init__(self, profile):
        """Clip the STL and set attributes on the C++
        GeometryGenerator object.
        """
        super().__init__()
        self.DebugPipelineDirectory = None
        self._profile = profile
        self.generator = self.genext.PolyDataGenerator()

    def _SetGeometryTypeProperties(self):
        self.generator.SetClippedSurface(self.ClippedSurface)

    def _Preprocess(self):
        # This will create the pipeline for the clipped surface
        clipper = Clipper(self._profile)
        # Note this is in coord system 2

        # Scale by the voxel size (to coord sys 3)
        scaleToLattice = vtkTransform()
        scale = 1.0 / self.VoxelSize
        scaleToLattice.Scale(scale, scale, scale)

        scaler = vtkTransformFilter()
        scaler.SetTransform(scaleToLattice)
        scaler.SetInputConnection(clipper.ClippedSurfaceSource.GetOutputPort())

        scaler.Update()

        originscaled, nSites = self._ComputeOriginScaled(scaler.GetOutput())
        self._profile.OriginMetres = Vector(originscaled * self.VoxelSizeMetres)
        self.SiteCounts = self._RoundUpSites(nSites)

        shiftOrigin = vtkTransform()
        shiftOrigin.Translate(-originscaled[0], -originscaled[1], -originscaled[2])
        shifter = vtkTransformFilter()
        shifter.SetTransform(shiftOrigin)
        shifter.SetInputConnection(scaler.GetOutputPort())

        shifter.Update()
        self.ClippedSurface = shifter.GetOutput()
        # If the attribute has been set, write all the stages of the
        # pipeline to the folder.
        if self.DebugPipelineDirectory is not None:
            sw = StageWriter(self.DebugPipelineDirectory)
            for alg in getpipeline(shifter):
                sw.WriteOutput(alg)

    @abstractmethod
    def _RoundUpSites(self, nSites):
        pass

    def _fixup_block(self):
        # Get the biggest of the sides
        nSites = int(nSites.max())
        # Now round this up to the next power of two
        nBits = nSites.bit_length()
        if nSites > 2 ** (nBits - 1):
            nSites = 2**nBits
            pass
        self.NumberOfLevels = nBits

        for i in range(3):
            # Now ensure this extra space is equally balanced before & after the
            # fluid region with the placement of the first site.
            bmin = SurfaceBoundsWorking[2 * i]
            bmax = SurfaceBoundsWorking[2 * i + 1]
            size = bmax - bmin
            extra = (nSites - 1) - size
            OriginWorking[i] = bmin - 0.5 * extra
            continue

        self.OriginWorking = OriginWorking
        self.CubeSize = nSites

    def _ComputeOriginScaled(self, surface):
        """Here we are working out the location of our domain's origin
        and the number of sites along each axis. This function assumes
        that we are already scaled to lattice units.  Sites will all
        have positions of:

        Origin + Index

        where: 0 <= Index[i] < nSites[i]

        We also require that there be at least one solid site outside the fluid
        sites. For the case of axis-aligned faces which are an integer number
        of VoxelSizes apart (e.g. synthetic datasets!) this can cause numerical
        issues for the classifier if all the points that are "outside" are very
        close to the surface so we further require that these sites are a
        little further from the bounding box of the PolyData.
        """
        SurfaceBounds = surface.GetBounds()
        Origin = np.zeros(3, dtype=float)
        nSites = np.zeros(3, dtype=np.uint)

        for i in range(3):
            # Bounds of the vtkPolyData
            bmin = SurfaceBounds[2 * i]
            bmax = SurfaceBounds[2 * i + 1]
            size = bmax - bmin

            nSites[i] = int(size)
            # Since int() truncates, we have:
            #         0 < size/VoxelSize - nSites < 1.
            # Hence we need nSites + 1 links and therefore nSites + 2 sites
            nSites[i] += 2

            # The extra distance from size to the distance from x[0] to x[nSites -1]
            extra = (nSites[i] - 1) - size

            # To avoid numerical problems with the classifier, ensure that the
            # sites just outside the fluid region are at least 1% of a VoxelSize
            # away.

            if extra < 0.01:
                # They weren't, so add one to the # sites and recalculate extra
                nSites[i] += 1
                extra = (nSites[i] - 1) - size
                pass

            # Now ensure this extra space is equally balanced before & after the
            # fluid region with the placement of the first site.
            Origin[i] = bmin - 0.5 * extra
            continue

        return Origin, nSites

    pass


class GmyPolyDataGenerator(GmyMixin, PolyDataGenerator):
    pass


class CylinderGenerator(GeometryGenerator, GmyMixin):
    def __init__(
        self,
        OutputGeometryFile,
        OutputXmlFile,
        VoxelSizeMetres,
        Axis,
        LengthMetres,
        RadiusMetres,
        InletPressure=None,
        OutletPressure=None,
    ):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        self.Axis = Axis
        self.LengthMetres = LengthMetres
        self.RadiusMetres = RadiusMetres
        self.InletPressure = InletPressure
        self.OutletPressure = OutletPressure

        self._profile = Profile()
        self._profile.StlFileUnitId = Profile._UnitChoices.index(metre)
        self._profile.VoxelSize = VoxelSizeMetres
        self._profile.OutputGeometryFile = OutputGeometryFile
        self._profile.OutputXmlFile = OutputXmlFile
        self._MakeIolets()

        self.generator = self.genext.CylinderGenerator()
        self._SetGeneratorProperties()

        self.generator.SetCylinderLength(LengthMetres / VoxelSizeMetres)
        self.generator.SetCylinderRadius(RadiusMetres / VoxelSizeMetres)
        self.generator.SetCylinderCentre(DoubleVector(0.0, 0.0, 0.0))
        self.generator.SetCylinderAxis(DoubleVector(*self.Axis))
        return

    def _MakeIolets(self):
        # Construct the Iolet structs
        inlet = Inlet()
        inlet.Centre = Vector(*(-0.5 * self.LengthMetres * n for n in self.Axis))
        inlet.Normal = Vector(*self.Axis)
        inlet.Radius = self.RadiusMetres
        if self.InletPressure is not None:
            inlet.Pressure = self.InletPressure
        self._profile.Iolets.append(inlet)

        outlet = Outlet()
        outlet.Centre = Vector(*(0.5 * self.LengthMetres * n for n in self.Axis))
        outlet.Normal = Vector(*(-n for n in self.Axis))
        outlet.Radius = self.RadiusMetres
        if self.OutletPressure is not None:
            outlet.Pressure = self.OutletPressure
        self._profile.Iolets.append(outlet)

        return

    pass


class SquareDuctGenerator(GeometryGenerator, GmyMixin):
    def __init__(
        self,
        OutputGeometryFile,
        OutputXmlFile,
        VoxelSizeMetres,
        OpenAxis,
        LengthVoxels,
        SideVoxels,
        InletPressure=None,
        OutletPressure=None,
    ):
        """Clip the STL and set attributes on the SWIG-proxied C++
        GeometryGenerator object.
        """
        GeometryGenerator.__init__(self)
        assert OpenAxis in (0, 1, 2)
        self.OpenAxis = OpenAxis
        self.LengthVoxels = LengthVoxels
        self.SideVoxels = SideVoxels
        self.Sizes = DoubleVector(SideVoxels, SideVoxels, SideVoxels)
        self.Sizes[OpenAxis] = LengthVoxels

        self.InletPressure = InletPressure
        self.OutletPressure = OutletPressure

        self._profile = Profile()
        self._profile.StlFileUnitId = Profile._UnitChoices.index(metre)
        self._profile.VoxelSize = VoxelSizeMetres
        self._profile.OutputGeometryFile = OutputGeometryFile
        self._profile.OutputXmlFile = OutputXmlFile
        self._MakeIolets()

        self.generator = self.genext.SquareDuctGenerator()
        self._SetGeneratorProperties()

        self.generator.SetOpenAxis(self.OpenAxis)
        lb = self.Sizes * -0.5
        self.generator.SetLowerBound(lb)
        ub = self.Sizes * 0.5
        self.generator.SetUpperBound(ub)
        return

    def _MakeIolets(self):
        # Construct the Iolet structs
        inlet = Inlet()
        c = [0.0, 0.0, 0.0]
        c[self.OpenAxis] = -0.5 * self.LengthVoxels * self._profile.VoxelSize
        inlet.Centre = Vector(*c)

        n = DoubleVector()
        n[self.OpenAxis] = 1.0
        inlet.Normal = Vector(n.x, n.y, n.z)

        inlet.Radius = self.SideVoxels * self._profile.VoxelSizeMetres
        if self.InletPressure is not None:
            inlet.Pressure = self.InletPressure
        self._profile.Iolets.append(inlet)

        outlet = Outlet()
        c = [0.0, 0.0, 0.0]
        c[self.OpenAxis] = 0.5 * self.LengthVoxels * self._profile.VoxelSize
        outlet.Centre = Vector(*c)

        n = DoubleVector()
        n[self.OpenAxis] = -1.0
        outlet.Normal = Vector(n.x, n.y, n.z)

        outlet.Radius = self.SideVoxels * self._profile.VoxelSizeMetres
        if self.OutletPressure is not None:
            outlet.Pressure = self.OutletPressure
        self._profile.Iolets.append(outlet)

        return

    pass


# TODO: organise this timer
class Timer:
    def __init__(self):
        self._running = False

    @staticmethod
    def get_hwm():
        # Only works on Linux
        if sys.platform != "linux":
            return None
        with open("/proc/{pid}/status".format(pid=os.getpid())) as stats:
            for line in stats:
                if line.startswith("VmHWM:"):
                    key, val, unit = line.split()
                    return (int(val), unit)

    def Start(self):
        assert not self._running
        self._running = True
        self._startTime = time.perf_counter()
        self._startMem = self.get_hwm()
        return

    def Stop(self):
        assert self._running
        self._stopTime = time.perf_counter()
        self._stopMem = self.get_hwm()
        self._running = False
        return

    def GetTime(self):
        if self._running:
            return time.perf_counter() - self._startTime
        return self._stopTime - self._startTime

    def GetMem(self):
        if self._startMem is None:
            return (math.nan, "B")
        stopMem = self.get_hwm() if self._running else self._stopMem
        assert self._startMem[1] == stopMem[1]
        return (stopMem[0] - self._startMem[0], stopMem[1])

    pass


class Clipper(object):
    """Clips the input STL file to the ROI and caps it."""

    def __init__(self, profile):
        """Create the generator for the supplied Model.Profile object.

        profile - the HlbGmyTool.Model.Profile object to use
        """
        # Copy the profile's _Args
        for k in profile._Args:
            setattr(self, k, getattr(profile, k))
            continue
        # Pull in the SurfaceSource too
        self.SurfaceSource = profile.StlReader

        self.ClippedSurfaceSource = self.ConstructClipPipeline()
        return

    def ConstructClipPipeline(self):
        """This constructs a VTK pipeline to clip the vtkPolyData read from
        the STL file against the Iolets. It also adds a scalar value to each
        polygon indicating which Iolet it represents, or -1 if it is just a
        wall, and computes polygon normals.

        Note that this does NOT EXECUTE the pipeline.
        """
        # seal any leaks first.
        closer = PolyDataCloser()
        closer.SetInputConnection(self.SurfaceSource.GetOutputPort())
        # Add the Iolet id -1 to all cells
        adder = IntegerAdder(Value=-1)
        adder.SetInputConnection(closer.GetOutputPort())
        # Have the name pdSource first point to the input, then loop
        # over IOlets, clipping and capping.
        pdSource = adder
        for i, iolet in enumerate(self.Iolets):
            capper = PolyDataClipCapAndLabeller(
                Value=i, Iolet=iolet, SeedPoint=self.SeedPoint
            )
            capper.SetInputConnection(pdSource.GetOutputPort())
            # Set the source of the next iteraction to the capped
            # surface producer.
            pdSource = capper
            continue

        # The following adds cell normals to the PolyData
        normer = vtkPolyDataNormals()
        normer.SetInputConnection(pdSource.GetOutputPort())
        normer.ComputeCellNormalsOn()
        normer.ComputePointNormalsOff()
        normer.SplittingOff()
        normer.ConsistencyOn()
        normer.AutoOrientNormalsOn()
        normer.NonManifoldTraversalOff()

        return normer

    pass


class IntegerAdder(vtkProgrammableFilter):
    """vtkFilter for adding an integer value to vtkPolyData's CellData. Use
    SetValue to set the value to be added, or supply the optional argument
    'Value' to the constructor.
    """

    def __init__(self, Value=-1):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        return

    def SetValue(self, val):
        self.Value = val
        return

    def GetValue(self):
        return self.Value

    def _Execute(self, *args):
        """Run the filter."""
        # Get the input
        input = self.GetPolyDataInput()
        values = vtkIntArray()
        values.SetNumberOfTuples(input.GetNumberOfCells())
        values.FillComponent(0, self.Value)

        # Get the out PD
        out = self.GetPolyDataOutput()
        out.CopyStructure(input)
        out.GetCellData().SetScalars(values)

        return

    pass


class PolyDataClipCapAndLabeller(vtkProgrammableFilter):
    """vtkFilter for clipping and capping a vtkPolyData surface, and labeling
    the cap with an integer cell data value.
    """

    def __init__(self, Value=None, Iolet=None, SeedPoint=None):
        self.SetExecuteMethod(self._Execute)
        self.Value = Value
        self.Iolet = Iolet
        self.SeedPoint = (SeedPoint.x, SeedPoint.y, SeedPoint.z)
        return

    def SetValue(self, val):
        self.Value = val
        return

    def GetValue(self):
        return self.Value

    def SetIolet(self, val):
        self.Iolet = val
        return

    def GetIolet(self):
        return self.Iolet

    def _Clip(self, pd):
        # The plane implicit function will be >0 for all the points in the positive side
        # of the plane (i.e. x s.t. n.(x-o)>0, where n is the plane normal and o is the
        # plane origin).
        plane = vtkPlane()
        plane.SetOrigin(self.Iolet.Centre.x, self.Iolet.Centre.y, self.Iolet.Centre.z)
        plane.SetNormal(self.Iolet.Normal.x, self.Iolet.Normal.y, self.Iolet.Normal.z)

        # The sphere implicit function will be >0 for all the points outside
        # the sphere.
        sphere = vtkSphere()
        sphere.SetCenter(self.Iolet.Centre.x, self.Iolet.Centre.y, self.Iolet.Centre.z)
        sphere.SetRadius(self.Iolet.Radius)

        # The VTK_INTERSECTION operator takes the maximum value of all the registered
        # implicit functions. This will result in the function evaluating to >0 for all
        # the points outside the sphere plus those inside the sphere in the positive
        # side of the plane.
        clippingFunction = vtkImplicitBoolean()
        clippingFunction.AddFunction(plane)
        clippingFunction.AddFunction(sphere)
        clippingFunction.SetOperationTypeToIntersection()

        clipper = vtkClipPolyData()
        clipper.SetInputData(pd)
        clipper.SetClipFunction(clippingFunction)

        # Filter to get part closest to seed point
        connectedRegionGetter = vtkPolyDataConnectivityFilter()
        connectedRegionGetter.SetExtractionModeToClosestPointRegion()
        connectedRegionGetter.SetClosestPoint(*self.SeedPoint)
        connectedRegionGetter.SetInputConnection(clipper.GetOutputPort())
        connectedRegionGetter.Update()
        return connectedRegionGetter.GetOutput()

    def _AddValue(self, pd):
        adder = IntegerAdder(Value=self.Value)
        adder.SetInputData(pd)
        adder.Update()
        return adder.GetOutput()

    def _Execute(self, *args):
        assert isinstance(self.Value, int)
        assert isinstance(self.Iolet, Iolet)
        input = self.GetPolyDataInput()
        clipped = self._Clip(input)

        clipped.BuildLinks(0)

        newPoints = vtkPoints()
        newPoints.DeepCopy(clipped.GetPoints())

        newPolys = vtkCellArray()
        newPolys.DeepCopy(clipped.GetPolys())

        newData = vtkIntArray()
        newData.DeepCopy(clipped.GetCellData().GetScalars())

        boundaryExtractor = vtkvmtkPolyDataBoundaryExtractor()
        boundaryExtractor.SetInputData(clipped)
        boundaryExtractor.Update()
        boundaries = boundaryExtractor.GetOutput()
        boundariesPointIdMap = boundaries.GetPointData().GetScalars()
        for i in range(boundaries.GetNumberOfCells()):
            boundary = vtkPolyLine.SafeDownCast(boundaries.GetCell(i))

            barycentre = [0.0, 0.0, 0.0]
            vtkvmtkBoundaryReferenceSystems.ComputeBoundaryBarycenter(
                boundary.GetPoints(), barycentre
            )

            barycentreId = newPoints.InsertNextPoint(barycentre)

            numberOfBoundaryPoints = boundary.GetNumberOfPoints()
            trianglePoints = vtkIdList()
            trianglePoints.SetNumberOfIds(3)

            for j in range(numberOfBoundaryPoints):
                trianglePoints.SetId(
                    0, boundariesPointIdMap.GetValue(boundary.GetPointId(j))
                )
                trianglePoints.SetId(1, barycentreId)
                trianglePoints.SetId(
                    2,
                    boundariesPointIdMap.GetValue(
                        boundary.GetPointId((j + 1) % numberOfBoundaryPoints)
                    ),
                )
                newPolys.InsertNextCell(trianglePoints)
                newData.InsertNextValue(self.Value)
                continue

            continue

        output = self.GetPolyDataOutput()
        output.SetPoints(newPoints)
        output.SetPolys(newPolys)
        output.GetCellData().SetScalars(newData)

        return

    pass


class PolyDataCloser(vtkProgrammableFilter):
    """Given an input vtkPolyData object, close any holes in the surface."""

    def __init__(self):
        self.SetExecuteMethod(self._Execute)

    def _Execute(self, *args):
        input = self.GetPolyDataInput()

        # Gets any edges of the mesh
        edger = vtkFeatureEdges()
        edger.BoundaryEdgesOn()
        edger.FeatureEdgesOff()
        edger.NonManifoldEdgesOff()
        edger.ManifoldEdgesOff()
        edger.SetInputData(input)

        # Converts the edges to a polyline
        stripper = vtkStripper()
        stripper.SetInputConnection(edger.GetOutputPort())
        stripper.Update()

        # Change the polylines into polygons
        boundaryPoly = vtkPolyData()
        boundaryPoly.SetPoints(stripper.GetOutput().GetPoints())
        boundaryPoly.SetPolys(stripper.GetOutput().GetLines())

        # Triangulate
        tri = vtkTriangleFilter()
        tri.SetInputData(boundaryPoly)
        tri.Update()

        # Join to the input
        merger = vtkAppendPolyData()
        merger.AddInputData(input)
        merger.AddInputData(tri.GetOutput())

        # Clean up by merging duplicate points
        cleaner = vtkCleanPolyData()
        cleaner.SetInputConnection(merger.GetOutputPort())
        cleaner.Update()

        # Set output
        output = self.GetPolyDataOutput()
        output.ShallowCopy(cleaner.GetOutput())
        return

    pass


# Debugging helpers


def getpreviousstage(algo):
    """Given a vtkAlgorithm, get the previous algorithm in its pipeline."""
    return algo.GetInputConnection(0, 0).GetProducer()


def getpipeline(last):
    """Return a list of all the algorithms in a pipeline, given the last one."""
    pipe = [last]
    while True:
        try:
            prev = getpreviousstage(pipe[0])
            pipe.insert(0, prev)
        except:
            break
        continue

    return pipe


class StageWriter(object):
    """For debugging, will easily let one write polydata."""

    def __init__(self, dir):
        if os.path.exists(dir):
            assert os.path.isdir(dir)
        else:
            os.mkdir(dir)
            pass
        self.dir = dir
        self.i = 0
        return

    def WriteOutput(self, stage, name=None):
        writer = vtkXMLPolyDataWriter()

        if isinstance(stage, vtkAlgorithm):
            # To ensure everything executes in order
            stage.Update()
            writer.SetInputConnection(stage.GetOutputPort())
        elif isinstance(stage, vtkPolyData):
            writer.SetInputData(stage)
        else:
            raise ValueError('Cannot cope with instances of "%s"' % type(stage))

        if name is None:
            fnPattern = "%02d.vtp"
        else:
            fnPattern = "%02d-" + name + ".vtp"
        filename = os.path.join(self.dir, fnPattern % self.i)
        self.i += 1
        writer.SetFileName(filename)
        writer.Write()
        return

    pass
