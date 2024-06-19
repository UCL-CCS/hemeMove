# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import h5py as h5
import numpy as np


class SectionTree(object):
    """Note all indexes are np.uint16"""

    def __init__(self, filename):
        self.filename = filename
        self.file = h5.File(filename, "r")
        self.root = self.file["hemelb/geometry"]

        assert self.root.attrs["version"] == 1, "Unknown version of OCT file!"

        self.levels = np.uint16(self.root.attrs["levels"])
        self.level_names = ["{:04d}".format(i) for i in range(self.levels + 1)]

        self.indices = [self.root[ln] for ln in self.level_names]

        return

    NA = ~np.uint64(0)
    one = np.uint16(1)

    @classmethod
    def local_offset(cls, i, j, k, lvl):
        """All args should be uint16"""
        xbit = (i >> lvl) & cls.one
        ybit = (j >> lvl) & cls.one
        zbit = (k >> lvl) & cls.one
        return (xbit << 2) | (ybit << 1) | zbit

    def find_index(self, i, j, k):
        """Args should be uint16

        return is uint64
        value of NA => does not exist
        """
        cur_level = self.levels
        # Root level (==nLevel) has implicit offset of zero
        cur_offset = np.uint64(0)
        while cur_level:
            # Get the local index
            lIndex = self.local_offset(i, j, k, cur_level - self.one)
            cur_offset = self.indices[cur_level][cur_offset + lIndex]
            if cur_offset == self.NA:
                break

            cur_level -= self.one
        return cur_offset

    def find_path(self, i, j, k, level=np.uint16(0)):
        cur_level = self.levels

        path = np.empty(self.levels + 1, dtype=np.uint64)
        # Root level (==nLevel) has implicit offset of zero
        # The rest are unknown
        path[:cur_level] = self.NA
        path[cur_level] = 0

        while cur_level > level:
            # Get the local index
            lIndex = self.local_offset(i, j, k, cur_level - self.one)

            new_off = path[cur_level - self.one] = self.indices[cur_level][
                path[cur_level] + lIndex
            ]
            if new_off == self.NA:
                break
            cur_level -= self.one
            continue
        return path

    def find_from(self, i, j, k, level, start_path, start_level):
        path = start_path.copy()
        cur_level = start_level
        while cur_level > level:
            # Get the local index
            lIndex = self.local_offset(i, j, k, cur_level - self.one)

            new_off = path[cur_level - self.one] = self.indices[cur_level][
                path[cur_level] + lIndex
            ]
            if new_off == self.NA:
                break
            cur_level -= self.one
            continue
        return path

    def get(self, section, i):
        """Get some data from the section at the requested index"""
        sg = self.root[section]
        count = sg["counts"][i]
        offset = sg["offsets"][i]
        return sg["data"][offset : (offset + count)]

    def link(self, i):
        return self.get("links", i)

    def wall(self, i):
        return self.get("wall_normals", i)

    pass
