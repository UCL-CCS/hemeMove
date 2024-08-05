# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os.path
import xdrlib
import numpy as np

from . import HemeLbMagicNumber

ExtractionMagicNumber = 0x78747204
MagicVersionLength = 12
TimeStepDataLength = 8


class FieldSpec:
    """Represent the data type of a single record in both XDR format and
    the native (fast) format of the machine.

    Fields are the point id, position in metres, index position and the fields
    specified in the extraction file.

    The specification of a single element is a five element tuple of:
        (field name,
         XDR data type,
         in-memory dtype,
         number of elements,
         offset from start of a "row",
         data offset,
         lattice->physical scale factor)

    """

    def __init__(self, memspec):
        # name, XDR dtype, in-memory dtype, length, offset
        self._filespec = [("grid", ">i4", np.uint32, (3,), 0, None, None)]

        self._memspec = memspec
        return

    def Append(self, name, length, xdrType, memType, d_off, scale):
        """Add a new field to the specification."""
        if length == 1:
            np_len = ()
        else:
            np_len = (length,)

        offset = self.GetRecordLength()
        self._filespec.append((name, xdrType, memType, np_len, offset, d_off, scale))
        return

    def GetMem(self):
        """Get the numpy datatype for the in-memory array."""
        return np.dtype(
            [
                (name, memType, length)
                for name, xdrType, memType, length, offset, d_off, scale in (
                    self._memspec + self._filespec
                )
            ]
        )

    def GetXdr(self):
        """Get the numpy datatype for the XDR file."""
        return np.dtype(
            [
                (name, xdrType, length)
                for name, xdrType, memType, length, offset, d_off, scale in self._filespec
            ]
        )

    def GetRecordLength(self):
        """Get the length of the record as stored in the XDR file."""
        return self.GetXdr().itemsize

    def __iter__(self):
        """Iterate over the file specification."""
        return iter(self._filespec)

    pass


class xtr_parser_base:
    def __init__(self, physical_units=True):
        if not physical_units:
            assert self.VERSION >= 6
        self.OutputPhysicalUnits = physical_units
        self.FieldHeaderStart = MagicVersionLength + self.MainHeaderLength

    def _load_field_header(self, xtr_file):
        xtr_file.seek(self.FieldHeaderStart)
        fieldHeader = xtr_file.read(self.FieldHeaderLength)
        assert (
            len(fieldHeader) == self.FieldHeaderLength
        ), "Did not read the correct length of the field header in extraction file '{}'".format(
            xtr_file.name
        )

        return xdrlib.Unpacker(fieldHeader)


class xtr_common_4_5:
    """Same main header for V4 & V5"""

    MainHeaderLength = 48

    def ParseMainHeader(self, xtr_file, xtr):
        # Read the correct number of bytes
        mainHeader = xtr_file.read(self.MainHeaderLength)
        assert (
            len(mainHeader) == self.MainHeaderLength
        ), "Did not read the correct length of the main header in extraction file '{}'".format(
            xtr_file.name
        )

        decoder = xdrlib.Unpacker(mainHeader)
        xtr.voxelSizeMetres = decoder.unpack_double()
        xtr.originMetres = np.array([decoder.unpack_double() for i in range(3)])

        self.siteCount = xtr.siteCount = decoder.unpack_uhyper()
        self.fieldCount = xtr.fieldCount = decoder.unpack_uint()
        self.FieldHeaderLength = decoder.unpack_uint()


class xtr_common_5_6:
    """Same field headers for V5 & 6"""

    TYPECODE_TYPE = [np.float32, np.float64, np.int32, np.uint32, np.int64, np.uint64]
    TYPECODE_STR = [">f4", ">f8", ">i4", ">u4", ">i8", ">u8"]
    UNPACK_TYPE = [
        lambda up: up.unpack_float,
        lambda up: up.unpack_double,
        lambda up: up.unpack_int,
        lambda up: up.unpack_uint,
        lambda up: up.unpack_hyper,
        lambda up: up.unpack_uhyper,
    ]

    def _parse_scale(self, decoder, tc):
        return None

    def ParseFieldHeader(self, xtr_file):
        decoder = self._load_field_header(xtr_file)

        self._fieldSpec = FieldSpec(
            [
                ("id", None, np.uint64, (), None, None, None),
                ("position", None, np.float32, (3,), None, None, None),
            ]
        )
        for iField in range(self.fieldCount):
            name = decoder.unpack_string().decode("ascii")
            length = decoder.unpack_uint()
            tc = decoder.unpack_uint()
            np_type = self.TYPECODE_TYPE[tc]
            n_offsets = decoder.unpack_uint()
            offsets = np.empty(n_offsets, dtype=np_type)
            for iOff in range(n_offsets):
                offsets[iOff] = self.UNPACK_TYPE[tc](decoder)()

            if n_offsets == 0:
                d_off = None
            elif n_offsets == 1:
                d_off = offsets[0]
            elif n_offsets == length:
                # Above condition will run for scalar fields
                d_off = offsets[np.newaxis, :]
            else:
                raise ValueError(
                    f"Invalid number of offsets in extraction file for field '{name}'"
                )
            scale = self._parse_scale(decoder, tc)
            self._fieldSpec.Append(
                name, length, self.TYPECODE_STR[tc], np_type, d_off, scale
            )

        return self._fieldSpec

    def ParseTimeStep(self, memoryMappedData):
        result = np.recarray(self.siteCount, dtype=self._fieldSpec.GetMem())

        for (
            name,
            xdrType,
            memType,
            length,
            offset,
            dataOffset,
            scale,
        ) in self._fieldSpec:
            filedata = memoryMappedData.getfield((xdrType, length), offset)
            memdata = getattr(result, name)
            memdata[:] = filedata[:]
            if dataOffset is not None:
                memdata += dataOffset
            if scale is not None:
                memdata *= scale
        return result


class ExtractedPropertyV4Parser(xtr_parser_base, xtr_common_4_5):
    """Only float32s, output in physical units"""

    VERSION = 4

    def ParseTimeStep(self, memoryMappedData):
        result = np.recarray(self.siteCount, dtype=self._fieldSpec.GetMem())

        for (
            name,
            xdrType,
            memType,
            length,
            offset,
            dataOffset,
            scale,
        ) in self._fieldSpec:
            data = memoryMappedData.getfield((xdrType, length), offset)
            setattr(result, name, self._recursiveAdd(data, dataOffset))
            continue
        return result

    def ParseFieldHeader(self, xtr_file):
        decoder = self._load_field_header(xtr_file)

        self._fieldSpec = FieldSpec(
            [
                ("id", None, np.uint64, (), None, None, None),
                ("position", None, np.float32, (3,), None, None, None),
            ]
        )
        self._dataOffset = [0]

        for iField in range(self.fieldCount):
            name = decoder.unpack_string().decode("ascii")
            length = decoder.unpack_uint()
            self._dataOffset.append(decoder.unpack_double())
            self._fieldSpec.Append(name, length, ">f4", np.float32)

        return self._fieldSpec

    def _recursiveAdd(self, data, operand):
        try:
            return [self._recursiveAdd(datum, operand) for datum in data]
        except TypeError:
            return data + operand
        pass


class ExtractedPropertyV5Parser(xtr_parser_base, xtr_common_4_5, xtr_common_5_6):
    """Variable types in physical units"""

    VERSION = 5
    pass


class ExtractedPropertyV6Parser(xtr_parser_base, xtr_common_5_6):
    """Variable types in lattice units"""

    VERSION = 6

    MainHeaderLength = 72

    def ParseMainHeader(self, xtr_file, xtr):
        # Read the correct number of bytes
        mainHeader = xtr_file.read(self.MainHeaderLength)
        assert (
            len(mainHeader) == self.MainHeaderLength
        ), "Did not read the correct length of the main header in extraction file '{}'".format(
            xtr_file.name
        )

        decoder = xdrlib.Unpacker(mainHeader)
        self.dx = xtr.voxelSizeMetres = decoder.unpack_double()
        self.dt = xtr.timeStepSeconds = decoder.unpack_double()
        self.dm = xtr.massScaleKg = decoder.unpack_double()
        self.origin_m = xtr.originMetres = np.array(
            [decoder.unpack_double() for i in range(3)]
        )
        self.ref_pressure_pa = decoder.unpack_double()
        self.siteCount = xtr.siteCount = decoder.unpack_uhyper()
        self.fieldCount = xtr.fieldCount = decoder.unpack_uint()
        self.FieldHeaderLength = decoder.unpack_uint()

    def _parse_scale(self, decoder, tc):
        scale = self.UNPACK_TYPE[tc](decoder)()
        if self.OutputPhysicalUnits:
            return None if scale == 0.0 else scale
        else:
            return None


class ExtractedProperty:
    """Represent the contents of a HemeLB property extraction file."""

    VersionHandlers = {
        4: ExtractedPropertyV4Parser,
        5: ExtractedPropertyV5Parser,
        6: ExtractedPropertyV6Parser,
    }

    def __init__(self, filename, physical_units=True):
        """Read the file's headers and determine how many times and which times
        have data available.
        """

        self.filename = filename
        self._file = open(filename, "rb")

        # Read the magic numbers and version
        magicEtc = self._file.read(MagicVersionLength)
        assert (
            len(magicEtc) == MagicVersionLength
        ), "Did not read the correct length of magic numbers in extraction file '{}'".format(
            self.filename
        )
        decoder = xdrlib.Unpacker(magicEtc)
        assert (
            decoder.unpack_uint() == HemeLbMagicNumber
        ), "Incorrect HemeLB magic number"
        assert (
            decoder.unpack_uint() == ExtractionMagicNumber
        ), "Incorrect extraction magic number"

        self.version = decoder.unpack_uint()
        try:
            self.ParserCls = self.VersionHandlers[self.version]
        except KeyError:
            raise RuntimeError("Incorrect extraction format version number")

        if not physical_units:
            # Unit conversion only supported for version >= 6
            assert self.version >= 6

        self.physical_units = physical_units
        self.parser = self.ParserCls()
        self.parser.ParseMainHeader(self._file, self)

        self._ReadFieldHeader()
        self._DetermineTimes()

        # At this point, we can close the file. All external access uses memory maps.
        self._file.close()
        return

    def _ReadFieldHeader(self):
        """Read the field headers. The main headers must have been read first.

        This also sets the _fieldSpec attribute, which is a list of tuples
        that give the name, datatype and shape of the record to be created.
        This is suitable to passing to numpy.dtype to create the recarray
        data type.
        """
        self._fieldSpec = self.parser.ParseFieldHeader(self._file)

        self._rowLength = self._fieldSpec.GetRecordLength()
        self._recordLength = TimeStepDataLength + self._rowLength * self.siteCount

        return

    def _DetermineTimes(self):
        """Examine the file to find out how many time steps worth of data and
        which times are contained within it.
        """
        filesize = os.path.getsize(self.filename)
        self._totalHeaderLength = (
            MagicVersionLength
            + self.parser.MainHeaderLength
            + self.parser.FieldHeaderLength
        )
        bodysize = filesize - self._totalHeaderLength
        assert bodysize % self._recordLength == 0, (
            "Extraction file appears to have partial record(s), residual %s / %s , bodysize %s"
            % (bodysize % self._recordLength, self._recordLength, bodysize)
        )
        nTimes = bodysize // self._recordLength

        times = np.zeros(nTimes, dtype=int)
        for iT in range(nTimes):
            pos = self._totalHeaderLength + iT * self._recordLength
            self._file.seek(pos)
            timeBuf = self._file.read(TimeStepDataLength)
            times[iT] = xdrlib.Unpacker(timeBuf).unpack_uhyper()
            continue

        assert np.all(
            np.argsort(times) == np.arange(len(times))
        ), "Times in extraction file are not monotonically increasing!"
        self.times = times

        return

    def GetByIndex(self, idx):
        """Get the fields by time index."""
        # Attempt to look up the index in the times array to catch any
        # IndexError that will be raised.
        t = self.times[idx]
        return self._LoadByIndex(idx)

    def GetByTimeStep(self, t):
        """Get the fields by time step."""
        idx = self.times.searchsorted(t)
        if self.times[idx] != t:
            raise IndexError(
                "Timestep {0} not in extraction file {1}".format(t, self.filename)
            )
        return self._LoadByIndex(idx)

    def GetFieldSpec(self):
        """Get the specification of all the fields we have"""
        return self._fieldSpec

    def _MemMap(self, idx):
        """Use numpy.memmap to make a single timestep's worth of data
        accessible through a numpy array.
        """
        # Figure out the start position of the data for this timestep in the
        # file. This is made up of
        #    - file headers
        #    - number of previous records * record length
        #    - stored timestep
        start = self._totalHeaderLength + idx * self._recordLength + TimeStepDataLength
        return np.memmap(
            self.filename,
            dtype=self._fieldSpec.GetXdr(),
            mode="r",
            offset=start,
            shape=(self.siteCount,),
        )

    def _LoadByIndex(self, idx):
        """Create a numpy record array with a single timestep of data.

        Fields are as specified in the file with the addition of
        """
        mapped = self._MemMap(idx)

        answer = self.parser.ParseTimeStep(mapped)

        answer.id = np.arange(self.siteCount)
        answer.position = self.voxelSizeMetres * answer.grid + self.originMetres
        return answer

    pass


if __name__ == "__main__":
    pass
