# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import xdrlib
import warnings

import numpy as np
from six.moves import range

from .. import HemeLbMagicNumber

SnapshotMagicNumber = 0x736E7004


def HemeLbSnapshot(filename):
    """Guess which file format we were given and use the correct class
    to open it.

    We have to handle a number of cases:

     - the original text format;

     - the XDR copy thereof, and

     - the updated (August 2011) version with format magic and version
       numbers and more metadata.
    """

    start = open(filename).read(8)
    reader = xdrlib.Unpacker(start)
    firstInt = reader.unpack_uint()

    if firstInt == HemeLbMagicNumber:
        assert reader.unpack_uint() == SnapshotMagicNumber
        cls = VersionedXdrSnapshot
    elif firstInt == 0 or firstInt == 1 or firstInt == 2:
        # It is the basic Xdr format that starts with the stablity flag
        cls = XdrSnapshotVersionOne
        # Maybe text? If so, the first character should be a '0', '1' or '2', followed by a newline
    elif (start[0] == "0" or start[0] == "1" or start == "2") and start[1] == "\n":
        cls = TextSnapshot
    else:
        raise ValueError('Cannot determine version of snapshot file "%s"' % filename)

    return cls(filename)


class BaseSnapshot(np.recarray):
    """Base class wrapping a HemeLB snapshot.

    Snap is basically a numpy record array with the following fields:

     - id (int) -- an id number (basically the index of the point in the
       file

     - position (3x float) -- the position in input space (m)

     - grid (3x int) -- the (x, y, z) coordinates in lattice units

     - pressure (float) -- the pressure in physical units (mmHg)

     - velocity (3x float) -- (x,y,z) components of the velocity field
       in physical units (m/s)

     - stress (float) -- the von Mises stress in physical units (Pa)

    It has a number of additional properties (see __readHeader for full details)

    """

    _raw_row = [
        ("id", int),
        ("position", float, (3,)),
        ("grid", int, (3,)),
        ("pressure", float),
        ("velocity", float, (3,)),
        ("stress", float),
    ]
    _readable_row = np.dtype(_raw_row[2:])
    row = np.dtype(_raw_row)

    _attrs = {
        "stable": None,
        "voxel_size": None,
        "origin": np.array([np.nan, np.nan, np.nan]),
        "bb_min": None,
        "bb_max": None,
        "bb_len": None,
        "voxel_count": None,
    }
    # header = len(_attrs)

    def __new__(cls, filename):
        """Create a new instance. Numpy array subclasses use this
        method instead of __init__ for initialization.

        """
        headerDict = cls._readHeader(filename)

        noindex = cls._load(filename, headerDict)
        index = np.recarray(shape=noindex.shape, dtype=cls.row)

        for el in cls._raw_row[2:]:
            key = el[0]
            index.__setattr__(key, noindex.__getattribute__(key))
            continue

        index.id = np.arange(len(noindex))
        try:
            index.position = cls._computePosition(index.grid, headerDict)
        except:
            index.position = np.nan
            pass

        obj = index.view(cls)

        # Set the attributes on the snapshot
        for headerField in headerDict:
            setattr(obj, headerField, headerDict[headerField])
            continue

        return obj

    def __array_finalize__(self, parent):
        """Numpy special method."""

        if parent is None:
            return
        for a in self._attrs:
            setattr(self, a, getattr(parent, a, self._attrs[a]))
            continue

        return

    pass


class PositionlessSnapshot(BaseSnapshot):
    """Base class for the original text snapshots and the XDR
    equivalent. These lack the data required to compute the positions
    of grid points. It is supplied through the coords.asc file
    generated by the old setuptool.
    """

    def computePosition(self, coordsFile):
        """Given the coordinate file from the segtool, calculate all
        the lattice positions' coordinates.

        """
        from os.path import exists

        if exists(coordsFile):
            from ...coordinates import Transformer

            trans = Transformer(coordsFile)

            self.position = 1e-3 * trans.siteToStl(self.grid + self.bb_min)
            return
        else:
            # The coords file is missing!
            warnings.warn(
                'Missing coordinates file "%s", assuming origin at [0,0,0]'
                % coordsFile,
                stacklevel=2,
            )
            self.position = (
                self.grid + self.bb_min
            ) * self.voxel_size  # + origin, but we'll just assume it's zero here.

    pass


class TextSnapshot(PositionlessSnapshot):
    """Read a text snapshot."""

    nHeaderLines = 6

    @classmethod
    def _readHeader(cls, filename):
        """Read the header lines, according to:
        0- Flag for simulation stability, 0 or 1
        1- Voxel size in physical units (units of m)
        2- vertex coords of the minimum bounding box with minimum values (x, y and z values)
        3- vertex coords of the minimum bounding box with maximum values (x, y and z values)
        4- #voxels within the minimum bounding box along the x, y, z axes (3 values)
        5- total number of fluid voxels

        """

        f = open(filename)
        stable = int(f.readline())
        voxel_size = float(f.readline())
        bb_min = np.array([int(x) for x in f.readline().split()])
        bb_max = np.array([int(x) for x in f.readline().split()])
        bb_len = np.array([int(x) for x in f.readline().split()])
        voxel_count = int(f.readline())

        return {
            "stable": stable,
            "voxel_size": voxel_size,
            "bb_min": bb_min,
            "bb_max": bb_max,
            "bb_len": bb_len,
            "voxel_count": voxel_count,
        }

    @classmethod
    def _load(cls, filename, header):
        return np.loadtxt(
            filename, skiprows=cls.nHeaderLines, dtype=cls._readable_row
        ).view(np.recarray)

    pass


class XdrVoxelFormatOneSnapshot(object):
    @classmethod
    def _load(cls, filename, header):
        # Skip past the header, slurp data, create XDR object
        f = open(filename)
        f.seek(cls._headerLengthBytes)
        reader = xdrlib.Unpacker(f.read())

        ans = np.recarray((header["voxel_count"],), dtype=cls._readable_row)
        # Read all the voxels.
        for i in range(header["voxel_count"]):
            ans[i] = (
                (reader.unpack_int(), reader.unpack_int(), reader.unpack_int()),
                reader.unpack_float(),
                (reader.unpack_float(), reader.unpack_float(), reader.unpack_float()),
                reader.unpack_float(),
            )

            continue
        reader.done()
        return ans

    pass


class XdrSnapshotVersionOne(PositionlessSnapshot, XdrVoxelFormatOneSnapshot):
    """Read an old-style XDR snapshot."""

    # int float 3x int 3x int 3x int int
    _headerLengthBytes = 4 + 8 + 3 * 4 + 3 * 4 + 3 * 4 + 4

    @classmethod
    def _readHeader(cls, filename):
        """Read the header lines, according to:
        0- Flag for simulation stability, 0 or 1
        1- Voxel size in physical units (units of m)
        2- vertex coords of the minimum bounding box with minimum values (x, y and z values)
        3- vertex coords of the minimum bounding box with maximum values (x, y and z values)
        4- #voxels within the minimum bounding box along the x, y, z axes (3 values)
        5- total number of fluid voxels

        """
        reader = xdrlib.Unpacker(open(filename).read(cls._headerLengthBytes))
        header = {}
        header["stable"] = reader.unpack_int()
        header["voxel_size"] = reader.unpack_double()

        header["bb_min"] = np.array(
            (reader.unpack_int(), reader.unpack_int(), reader.unpack_int())
        )
        header["bb_max"] = np.array(
            (reader.unpack_int(), reader.unpack_int(), reader.unpack_int())
        )
        header["bb_len"] = np.array(
            (reader.unpack_int(), reader.unpack_int(), reader.unpack_int())
        )
        header["voxel_count"] = reader.unpack_int()
        return header

    pass


class XdrSnapshotVersionTwo(BaseSnapshot, XdrVoxelFormatOneSnapshot):
    """Read snapshots for the updated format as for August 2011."""

    _headerLengthBytes = 80
    VersionNumber = 2

    @classmethod
    def _readHeader(cls, filename):
        """Read the header lines, according to description in Code/io/formats/snapshot.h"""
        reader = xdrlib.Unpacker(open(filename).read(cls._headerLengthBytes))
        header = {}
        assert reader.unpack_uint() == HemeLbMagicNumber
        assert reader.unpack_uint() == SnapshotMagicNumber
        assert reader.unpack_uint() == cls.VersionNumber
        bodyStart = reader.unpack_uint()
        assert bodyStart == cls._headerLengthBytes

        header["stable"] = reader.unpack_int()
        header["voxel_size"] = reader.unpack_double()
        header["origin"] = np.array(
            (reader.unpack_double(), reader.unpack_double(), reader.unpack_double())
        )

        header["bb_min"] = np.array(
            (reader.unpack_int(), reader.unpack_int(), reader.unpack_int())
        )
        header["bb_max"] = np.array(
            (reader.unpack_int(), reader.unpack_int(), reader.unpack_int())
        )
        header["bb_len"] = header["bb_max"] - header["bb_min"] + 1

        header["voxel_count"] = reader.unpack_int()
        return header

    @classmethod
    def _computePosition(cls, grid, header):
        return (grid + header["bb_min"]) * header["voxel_size"] + header["origin"]

    pass


def VersionedXdrSnapshot(filename):
    """Examine the file and dispatch to the appropriate constructor."""
    # Need the two magic numbers and the version number, i.e. 12 bytes
    reader = xdrlib.Unpacker(open(filename).read(12))

    assert reader.unpack_uint() == HemeLbMagicNumber
    assert reader.unpack_uint() == SnapshotMagicNumber
    version = reader.unpack_uint()
    if version == 2:
        return XdrSnapshotVersionTwo(filename)

    raise ValueError('Unknown version number (%d) in file "%s"' % (version, filename))
