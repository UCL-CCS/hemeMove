#! /usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import sys
from ..parsers.extraction import ExtractedProperty


def unpack(filename, stream=sys.stdout):
    propFile = ExtractedProperty(filename)

    print('# Dump of file "{}"'.format(filename), file=stream)
    print("# File has {} sites.".format(propFile.siteCount), file=stream)
    print("# File has {} fields:".format(propFile.fieldCount), file=stream)
    for name, xdrType, memType, length, offset, d_off, scale in propFile._fieldSpec:
        print('#     "{0}", length {1}'.format(name, length), file=stream)
    print("# Geometry origin = {} m".format(propFile.originMetres), file=stream)
    print("# Voxel size = {} m".format(propFile.voxelSizeMetres), file=stream)

    header = "# " + ", ".join(
        name
        for name, xdrType, memType, length, offset, d_off, scale in propFile._fieldSpec
    )
    print(header, file=stream)

    for t in propFile.times:
        fields = propFile.GetByTimeStep(t)
        print("# Timestep {:d}".format(t), file=stream)

        for row in fields:
            print(
                ", ".join(
                    str(getattr(row, name))
                    for name, xdrType, memType, length, offset, d_off, scale in propFile._fieldSpec
                ),
                file=stream,
            )

        print("", file=stream)


def main():
    unpack(sys.argv[1])
