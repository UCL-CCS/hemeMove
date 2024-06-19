# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import pytest

import filecmp
import os.path
import shutil
from xml.etree import ElementTree

import numpy as np
from numpy.linalg import norm

from HlbGmyTool.Model import OutputGeneration
from HlbGmyTool.Model.Profile import Profile
from HlbGmyTool.Model.Vector import Vector
from HlbGmyTool.Model.Iolets import Iolet
from hlb.converters.Oct2Gmy import oct2gmy
from hlb.parsers.geometry.simple import ConfigLoader
from hlb.parsers.geometry.self_consistency import CheckingLoader
from hlb.parsers.geometry.simple import ConfigLoader
from hlb.parsers.geometry import generic
from hlb.parsers.octree import SectionTree
from hlb.utils.xml_compare import XmlChecker
import fixtures

dataDir = os.path.join(os.path.split(__file__)[0], "data")


def vec2np(v):
    ans = np.empty(3, dtype=float)
    ans[0] = v.x
    ans[1] = v.y
    ans[2] = v.z
    return ans


class TestPolyDataGenerator:
    def test_regression(self, tmpdir):
        """Generate a gmy from a stored profile and check that the output is
        identical.
        """
        proFileName = os.path.join(dataDir, "test.pr2")

        p = Profile()
        p.LoadFromFile(proFileName)
        # Change the output to the tmpdir
        basename = tmpdir.join("test").strpath
        outGmyFileName = basename + ".gmy"
        outXmlFileName = basename + ".xml"
        p.OutputGeometryFile = outGmyFileName
        p.OutputXmlFile = outXmlFileName

        generator = OutputGeneration.GmyPolyDataGenerator(p)
        generator.Execute()

        ldr = CheckingLoader(outGmyFileName)
        ldr.Load()

        ref_ldr = ConfigLoader(os.path.join(dataDir, "test.gmy"))
        ref_ldr.Load()
        ref_dom = ref_ldr.Domain

        test_ldr = ConfigLoader(outGmyFileName)
        test_ldr.Load()
        test_dom = test_ldr.Domain

        # Per-block data length must be identical
        assert np.all(
            ref_ldr.BlockUncompressedDataLength == test_ldr.BlockUncompressedDataLength
        )
        # Same for fluid site counts
        assert np.all(ref_dom.BlockFluidSiteCounts == test_dom.BlockFluidSiteCounts)
        nblocks = len(ref_dom.BlockFluidSiteCounts)

        for i_b in range(nblocks):
            # Blocks should be the same
            ref_blk = ref_dom.Blocks[i_b]
            test_blk = test_dom.Blocks[i_b]
            assert type(ref_blk) == type(ref_blk)
            if isinstance(ref_blk, generic.AllSolidBlock):
                continue
            for i_s in range(ref_dom.BlockSize**3):
                rs = ref_blk.Sites[i_s]
                ts = test_blk.Sites[i_s]
                # Sites must match
                assert rs.Type == ts.Type
                if rs.IntersectionType is None:
                    assert ts.IntersectionType is None
                else:
                    assert np.all(rs.IntersectionType == ts.IntersectionType)

                if rs.IntersectionDistance is None:
                    assert ts.IntersectionDistance is None
                else:
                    assert np.allclose(
                        rs.IntersectionDistance,
                        ts.IntersectionDistance,
                        rtol=0,
                        atol=5e-5,
                    )

                assert np.all(rs.IOletIndex == ts.IOletIndex)
                assert rs.WallNormalAvailable == ts.WallNormalAvailable
                if rs.WallNormalAvailable:
                    assert np.allclose(rs.WallNormal, ts.WallNormal, atol=1e-6)

        # XML output matches also
        xmlChecker = XmlChecker.from_path(os.path.join(dataDir, "test.xml"))
        xmlChecker.check_path(outXmlFileName)

    def test_cube(self, tmpdir):
        """Generate a gmy from a simple cubic profile and check the output"""
        cube = fixtures.cube(tmpdir)
        cube.StlFileUnitId = 0
        cube.VoxelSize = 0.11
        generator = OutputGeneration.GmyPolyDataGenerator(cube)
        generator.Execute()
        # Load back the resulting geometry file and assert things are as
        # expected
        checker = CubeTestingGmyParser(cube.OutputXmlFile, cube.VoxelSize)
        checker.Load()

        fluid_sites = checker.Domain.BlockFluidSiteCounts.sum()
        block_count = len(checker.Domain.Blocks)
        block_size = checker.Domain.BlockSize
        sites = block_count * block_size**3
        assert sites == 4096
        assert fluid_sites == 729
        assert sites != fluid_sites
        # # Now, turn on the skip-non-intersecting-blocks optimisation, and
        # # assert same result
        # generator.skipNonIntersectingBlocks = True
        # generator.Execute()
        # checker_skip_nonintersecting = CubeTestingGmyParser(
        #     cube.OutputXmlFile, cube.VoxelSize)
        # checker_skip_nonintersecting.Load()
        # fluid_sites_nonintersecting = sum(
        #     checker_skip_nonintersecting.Domain.BlockFluidSiteCounts)
        # assert(fluid_sites_nonintersecting == fluid_sites)

    def test_cube_normals(self, tmpdir):
        """Generate a gmy from a simple cubic profile and check the computed
        normals.
        """
        cube = fixtures.cube(tmpdir)
        cube.VoxelSize = 0.23
        cube.StlFileUnitId = 0

        """The default VTK cube has 1m edges and it is centred at the origin 
        of coordinates. We place the inlet and the outlet at the faces 
        perpendicular to the z axis.
        """
        inlet = Iolet(
            Name="inlet",
            Centre=Vector(0.0, 0.0, -0.5),
            Normal=Vector(0.0, 0.0, -1.0),
            Radius=np.sqrt(2) / 2,
        )
        outlet = Iolet(
            Name="outlet",
            Centre=Vector(0.0, 0.0, 0.5),
            Normal=Vector(0.0, 0.0, 1.0),
            Radius=np.sqrt(2) / 2,
        )
        cube.Iolets = [inlet, outlet]

        generator = OutputGeneration.GmyPolyDataGenerator(cube)
        # generator.skipNonIntersectingBlocks = True
        generator.Execute()

        """Load back the resulting geometry file and assert things are as 
        expected
        """
        checker = CubeNormalsTestingGmyParser(cube.OutputXmlFile, cube.VoxelSize)
        checker.Load()

    def test_cylinder(self, tmpdir):
        """Generate a gmy from a cylinder and check the output.

        The cylinder is 1 length unit long, radius 0.5, aligned with the
        y-axis, and centred at the origin of coordinates.
        """
        dx = 0.11
        cylinder = fixtures.cylinder(tmpdir)
        cylinder.VoxelSize = dx
        cylinder.StlFileUnitId = 0

        inlet = Iolet(
            Name="inlet",
            Centre=Vector(0.0, -0.45, 0.0),
            Normal=Vector(0.0, -1.0, 0.0),
            Radius=1.0,
        )
        outlet = Iolet(
            Name="outlet",
            Centre=Vector(0.0, 0.45, 0.0),
            Normal=Vector(0.0, 1.0, 0.0),
            Radius=1.0,
        )
        cylinder.Iolets = [inlet, outlet]
        L = 0.9
        R = 0.5

        generator = OutputGeneration.GmyPolyDataGenerator(cylinder)

        generator.Execute()
        # Load back the resulting geometry file and assert things are as
        # expected
        checker = CylinderTestingGmyParser(
            cylinder.OutputXmlFile,
            cylinder.VoxelSize,
            np.array([0.0, 1.0, 0.0]),
            L,
            R,
        )
        checker.Load()

        fluid_sites = checker.Domain.BlockFluidSiteCounts.sum()
        block_count = len(checker.Domain.Blocks)
        block_size = checker.Domain.BlockSize
        sites = block_count * block_size**3
        # assert(sites==4096)
        assert fluid_sites == 621
        assert sites != fluid_sites
        # # Now, turn on the skip-non-intersecting-blocks optimisation, and
        # # assert same result
        # generator.skipNonIntersectingBlocks = True
        # generator.Execute()
        # checker_skip_nonintersecting = CylinderTestingGmyParser(
        #     cylinder.OutputGeometryFile, cylinder.VoxelSize,
        #     np.array([0.0, 1.0, 0.0]), 1.0, 0.5)
        # checker_skip_nonintersecting.Load()
        # fluid_sites_nonintersecting = sum(
        #     checker_skip_nonintersecting.Domain.BlockFluidSiteCounts)
        # assert(fluid_sites_nonintersecting == fluid_sites)


class TestCylinderGenerator:
    @pytest.mark.parametrize(("randSeed",), [(828,), (341,), (1432,)])
    def test_regression(self, tmpdir, randSeed):
        """Generate a small cylinder GMY with a random orientation. Then check
        that the output is correct by running it through a custom subclass of
        a ConfigLoader.
        """

        basename = tmpdir.join("cyl")
        OutputGeometryFile = basename.strpath + ".gmy"
        OutputXmlFile = basename.strpath + ".xml"
        VoxelSizeMetres = 0.1

        rng = np.random.RandomState(randSeed)
        Axis = rng.normal(size=(3,))
        Axis /= np.sqrt(np.dot(Axis, Axis))

        LengthMetres = 1.32
        RadiusMetres = 0.74

        # generator = OutputGeneration.CylinderGenerator(
        #     OutputGeometryFile, OutputXmlFile,
        #     VoxelSizeMetres, Axis, LengthMetres,
        #     RadiusMetres)
        # generator.Execute()

        # checker = CylinderTestingGmyParser(OutputGeometryFile, VoxelSizeMetres,
        #                                    Axis, LengthMetres, RadiusMetres)
        # checker.Load()
        return

    pass


class BaseTestingGmyParser(ConfigLoader):
    def __init__(self, filename, VoxelSize):
        ConfigLoader.__init__(self, filename)
        self.VoxelSize = VoxelSize

    def OnEndHeader(self):
        assert self.Domain.VoxelSize == self.VoxelSize

    def OnEndSite(self, block, site):
        if site.IsFluid:
            assert self.IsInside(site.Position)
        else:
            assert not self.IsInside(site.Position)
        return

    def OnEndBlock(self, bIdx, bIjk):
        self.Domain.DeleteBlock(bIdx)
        return

    pass


class CylinderTestingGmyParser(BaseTestingGmyParser):
    def __init__(self, filename, VoxelSize, Axis, Length, Radius):
        BaseTestingGmyParser.__init__(self, filename, VoxelSize)
        self.Axis = Axis
        self.Length = Length
        self.Radius = Radius
        return

    def IsInside(self, x):
        xDOTn = np.dot(x, self.Axis)
        if xDOTn < -0.5 * self.Length or xDOTn > 0.5 * self.Length:
            return False
        perp = x - xDOTn * self.Axis
        if np.dot(perp, perp) > self.Radius**2:
            return False
        return True

    def OnEndSite(self, block, site):
        BaseTestingGmyParser.OnEndSite(self, block, site)
        assert site.IsEdge == site.WallNormalAvailable
        is_cylinder_end = np.any(
            site.IntersectionType == site.INLET_INTERSECTION
        ) or np.any(site.IntersectionType == site.OUTLET_INTERSECTION)
        if site.IsEdge and not is_cylinder_end:
            assert np.any(site.IntersectionType == site.WALL_INTERSECTION)
            axis_perpendicular_at_site = (
                site.Position - np.dot(site.Position, self.Axis) * self.Axis
            )
            axis_perpendicular_at_site /= norm(axis_perpendicular_at_site)
            """ (ticket #597) 0.02 is an arbitrary tolerance for how accurate 
            the wall normal estimates are.
            """
            assert (
                np.absolute(1 - np.dot(site.WallNormal, axis_perpendicular_at_site))
                < 0.02
            )


class CubeTestingGmyParser(BaseTestingGmyParser):
    def IsInside(self, position):
        result = all(component < 0.5 and component > -0.5 for component in position)
        return result


class CubeNormalsTestingGmyParser(CubeTestingGmyParser):

    ValidNormals = np.array([(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)], dtype=float)

    def OnEndSite(self, block, site):
        CubeTestingGmyParser.OnEndSite(self, block, site)
        assert site.IsEdge == site.WallNormalAvailable
        is_cube_edge = (site.Index[0] in [1, 4]) and (site.Index[1] in [1, 4])
        if site.IsEdge and not is_cube_edge:
            assert np.any(site.IntersectionType == site.WALL_INTERSECTION)
            assert np.any(np.all(self.ValidNormals == site.WallNormal, axis=1))
