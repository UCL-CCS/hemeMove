# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
"""Define a filter to add the data in a HemeLB extracted property file to the geometry
defined by a HemeLB geometry file.

This module can be run as a script on the command line to convert an extraction to
the corresponding .vtu (VTk Unstructured grid XML file), for usage, run the
script with no arguments. 
"""
from xml.etree import ElementTree as et

import numpy as np
import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.util import numpy_support

from ..utils import MatchCorresponding
from ..parsers.extraction import ExtractedProperty


class ExtractedPropertyUnstructuredGridReader(VTKPythonAlgorithmBase):
    """VTK-style filter for reading HemeLB extracted property files as VTK
    data into the geometry provided as input. This input must be that derived
    from the HemeLB geometry file used by the simulation from which the
    extraction file comes and be as output by GmyUnstructuredGridReader (e.g.
    no scaling), as this class uses the positions to match points.

    The vtkUnstructuredGrid will have that subset of the points and cells from
    the input that contain data. Cell data will be as in the extraction file.
    The fields take the names and units given in the extracted property file.
    Position units are metres. The object has no point data.
    """

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self,
            nInputPorts=1,
            inputType="vtkUnstructuredGrid",
            nOutputPorts=1,
            outputType="vtkUnstructuredGrid",
        )
        self.Extracted = None
        self.Time = None
        self.Skeleton = None
        self.OutputCellIdsByInputIndex = None
        self.InputMTime = -1
        return

    def SetExtraction(self, extracted):
        """The parsed extracted property file."""
        self.Extracted = extracted
        self.Skeleton = None
        self.OutputCellIdsByInputIndex = None
        self.Modified()
        return

    def SetTime(self, time):
        """The timestamp to read properties for."""
        time = int(time)
        if self.Time != time:
            self.Time = time
            self.Modified()
        return

    def _CreateSkeleton(self, input_ug):
        """Create the structure of the output vtkUnstructuredGrid and a map
        from the index of a point in the extraction file to the corresponding
        cellId in the skeleton.

        This method should only be called if the extraction object this
        instance is working on has changed since the last time this method was
        called.
        """
        # Get the centres as these should match the extracted property positions
        centers = vtk.vtkCellCenters()
        centers.SetInputData(input_ug)
        centers.Update()
        # Use this to find the cell ID for each point.
        locator = vtk.vtkOctreePointLocator()
        locator.SetDataSet(centers.GetOutput())
        locator.BuildLocator()
        # Should be fine enough
        locator.SetTolerance(0.1 * self.Extracted.voxelSizeMetres)

        # Get the first set of extracted data from the file. We don't care
        # which as we only want to use the positions.
        extracted_positions = self.Extracted.GetByIndex(0).position
        nExtractedPoints = len(extracted_positions)

        # Make a list of the cell ids to keep; i.e. the cell in the input
        # (whole geometry) with ID cellIdsGmy[i] contains  the point given by
        # extracted_positions[i]
        gmyCellIdsByInput = vtk.vtkIdTypeArray()
        gmyCellIdsByInput.SetNumberOfComponents(1)
        gmyCellIdsByInput.SetNumberOfTuples(nExtractedPoints)

        gmyCellIdsByInputNp = numpy_support.vtk_to_numpy(gmyCellIdsByInput)

        for i, point in enumerate(extracted_positions):
            # Get cell in geometry corresponding to the point
            cellId = locator.FindClosestPoint(point)

            if cellId == -1:
                raise ValueError("Can't find cell for point at " + str(point))

            gmyCellIdsByInputNp[i] = cellId

        # Make an object to select only the cell ids we want
        selector = vtk.vtkSelectionNode()
        selector.SetFieldType(vtk.vtkSelectionNode.CELL)
        selector.SetContentType(vtk.vtkSelectionNode.INDICES)
        selector.SetSelectionList(gmyCellIdsByInput)

        # Make an object to hold the selector
        selectors = vtk.vtkSelection()
        selectors.AddNode(selector)

        # Perform the selection
        extractSelection = vtk.vtkExtractSelection()
        extractSelection.SetInputData(0, input_ug)
        extractSelection.SetInputData(1, selectors)
        extractSelection.Update()

        gmyCellIdsByOutput = (
            extractSelection.GetOutput().GetCellData().GetArray("vtkOriginalCellIds")
        )

        gmyCellIdsByOutputNp = numpy_support.vtk_to_numpy(gmyCellIdsByOutput)

        self.OutputCellIdsByInputIndex = MatchCorresponding(
            gmyCellIdsByOutputNp, gmyCellIdsByInputNp
        )

        # Create the skeleton data object and only copy in the structure.
        self.Skeleton = vtk.vtkUnstructuredGrid()
        self.Skeleton.CopyStructure(extractSelection.GetOutput())
        return

    def RequestData(self, request, inInfo, outInfo):
        """Method that actually does the reading. Called by the VTK
        API.
        """
        inp = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        # Check if we have a skeleton and if so, if it is up-to-date.
        if self.Skeleton is None or (inp.GetMTime() > self.InputMTime):
            self._CreateSkeleton(inp)
            self.InputMTime = inp.GetMTime()

        # Copy the structure to output.
        output = vtk.vtkUnstructuredGrid.GetData(outInfo)
        output.CopyStructure(self.Skeleton)

        # Get the data we want to visualise, Note that this is likely to
        # involve a subset of the points in the geometry
        extracted_data = self.Extracted.GetByTimeStep(self.Time)
        nExtractedPoints = len(extracted_data)

        # Create the arrays to store field data
        field_dict = {}

        for (
            name,
            xdrType,
            memType,
            length,
            offset,
            d_off,
            scale,
        ) in self.Extracted.GetFieldSpec():
            # Skip the grid (i.e. 3d-index)
            if name == "grid":
                continue

            # Set up the VTK array
            field = vtk.vtkDoubleArray()
            length = length[0] if len(length) else 1
            field.SetNumberOfComponents(length)
            field.SetNumberOfTuples(nExtractedPoints)
            field.SetName(name)

            # Insert it into the dictionary
            field_dict[(name, length)] = field

        # Copy the data into the correct position in the output array and add
        # the array to the output.
        for (field_name, field_length), field in field_dict.items():
            # fieldArray is a view into the data stored in field
            fieldArray = numpy_support.vtk_to_numpy(field)
            data = getattr(extracted_data, field_name)

            if field_length == 6:
                # VTK assumes the order XX, YY, ZZ, XY, YZ, XZ for symmetric
                # tensors stored in compressed format. However HemeLB outputs
                # XX XY XZ YY YZ ZZ. Reorder elements.
                fieldArray[self.OutputCellIdsByInputIndex] = data[:, [0, 3, 5, 1, 4, 2]]
            else:
                fieldArray[self.OutputCellIdsByInputIndex] = data

            output.GetCellData().AddArray(field)

        return 1


def WritePVDFile(timestepToFileMap, baseFilename):
    # Root element defining metadata
    vtkFile = et.Element(
        "VTKFile",
        {
            "type": "Collection",
            "version": "0.1",
            "byte_order": "LittleEndian",
            "compressor": "vtkZLibDataCompressor",
        },
    )
    # Collection element containing individual time steps
    collection = et.SubElement(vtkFile, "Collection")
    # Write out each time step and associated file
    for timestep, filename in timestepToFileMap.items():
        et.SubElement(
            collection,
            "DataSet",
            {"timestep": str(timestep), "group": "", "part": "0", "file": filename},
        )

    # Create ElementTree object and write it to disk
    xmlTree = et.ElementTree(vtkFile)
    xmlTree.write("{}.pvd".format(baseFilename))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create a VTK Unstructured grid with geometry of the "
        "first argument and the data in the extracted property file in all later arguments."
    )
    parser.add_argument(
        "geometry",
        nargs=1,
        help="the geometry, either a HemeLB .gmy file "
        "or a derived VTK unstructured grid",
    )
    parser.add_argument(
        "extracted",
        nargs=argparse.ONE_OR_MORE,
        help="extracted property file to convert to VTK, output will"
        ' be put in the input basename + ".vtu"',
    )

    args = parser.parse_args()
    geometry = args.geometry[0]
    import os.path

    base, ext = os.path.splitext(geometry)
    if ext == ".gmy":
        from .GmyUnstructuredGridReader import GmyUnstructuredGridReader

        reader = GmyUnstructuredGridReader()
    elif ext == ".vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        print('Cannot infer reader time from file extension "%s"' % ext)
        parser.print_help()
        raise SystemExit()

    reader.SetFileName(geometry)

    # For each extracted file given, parse the file...
    for extractedFile in args.extracted:
        extraction = ExtractedProperty(extractedFile)
        converter = ExtractedPropertyUnstructuredGridReader()
        converter.SetInputConnection(reader.GetOutputPort())
        converter.SetExtraction(extraction)

        # Create a writer...
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputConnection(converter.GetOutputPort())

        base, ext = os.path.splitext(extractedFile)

        # For each timestamp in that file
        timestepToFileMap = {}
        for time in extraction.times:
            # Choose the time
            converter.SetTime(time)

            # And give the output a unique filename based on the extracted
            # property file name and the timestamp.
            filename = "%s_%s.vtu" % (base, str(time))
            writer.SetFileName(filename)
            # trigger the update
            writer.Write()

            # Keep track of the files written for each time step
            timestepToFileMap[time] = filename

        # Write out file combining all time steps for easier visualisation
        WritePVDFile(timestepToFileMap, base)
