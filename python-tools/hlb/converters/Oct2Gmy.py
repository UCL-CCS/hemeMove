# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import argparse
import os.path
import xdrlib
import zlib

import numpy as np

from ..parsers.octree import SectionTree
from ..parsers.geometry import HemeLbMagicNumber, GeometryMagicNumber


class GmyBuilder(object):
    GMY_VERSION = 4
    BLOCK_SIZE = 8

    def __init__(self, fn, size):
        self.filename = fn

        self.size_sites = np.array((size, size, size), dtype=np.uint)
        self.size_blocks = self.size_sites // self.BLOCK_SIZE
        self.n_blocks = np.prod(self.size_blocks)

    def __enter__(self):
        """Context manager style access.

        This writes the preamble and blank headers.

        Call the returned object to add a block

        Exit method will clean up and rewrite correct headers.
        """
        self.outfile = open(self.filename, "wb")
        self._write_preamble()
        self.header_start_pos = self.outfile.tell()
        self.header = np.zeros((self.n_blocks, 3), dtype=np.uint)
        self._write_header()

        return self

    def __exit__(self, *args):
        """Rewrite the header and close file"""
        self._write_header()
        self.outfile.close()
        return

    def __iter__(self):
        self.current_block = np.zeros(3, dtype=np.uint16)
        # Note that Python immediately calls next on an iterator
        self.current_block[2] -= 1
        self.current_i = -1
        return self

    def __next__(self):
        """Next block"""
        self.current_i += 1
        self.current_block[2] += 1
        if self.current_block[2] == self.size_blocks[2]:
            self.current_block[1] += 1
            if self.current_block[1] == self.size_blocks[1]:
                self.current_block[0] += 1
                if self.current_block[0] == self.size_blocks[0]:
                    raise StopIteration
                self.current_block[1] = 0
            self.current_block[2] = 0

        return BlockBuilder(self, self.current_block)

    def _write_block(self, nsites, uncompressed):
        if nsites:
            compressed = zlib.compress(uncompressed)
            self.header[self.current_i] = nsites, len(compressed), len(uncompressed)
            self.outfile.write(compressed)

    def _write_preamble(self):
        p = xdrlib.Packer()

        p.pack_uint(HemeLbMagicNumber)
        p.pack_uint(GeometryMagicNumber)
        p.pack_uint(self.GMY_VERSION)

        for x in self.size_blocks:
            p.pack_uint(x)
        p.pack_uint(self.BLOCK_SIZE)

        p.pack_uint(0)  # padding

        self.outfile.write(p.get_buffer())
        return

    def _write_header(self):
        p = xdrlib.Packer()
        for block in self.header:
            p.pack_uint(block[0])
            p.pack_uint(block[1])
            p.pack_uint(block[2])

        self.outfile.seek(self.header_start_pos)
        self.outfile.write(p.get_buffer())
        return

    pass


class BlockBuilder(object):
    def __init__(self, gb, idx):
        self.gmy_builder = gb
        self.block_index = idx.copy()
        self.site_zero_index = self.block_index << 3
        return

    def __enter__(self):
        self.site_count = 0
        self.packer = xdrlib.Packer()

        return self

    def __exit__(self, *args):
        self.gmy_builder._write_block(self.site_count, self.packer.get_buffer())
        pass

    pass


def oct2gmy(octfilename, gmyfilename=None):
    if gmyfilename is None:
        base, ext = os.path.splitext(octfilename)
        gmyfilename = base + ".gmy"

    tree = SectionTree(octfilename)
    NA = SectionTree.NA
    box_size = 2**tree.levels

    with GmyBuilder(gmyfilename, box_size) as gb:
        for bb in gb:
            with bb:
                s0 = bb.site_zero_index
                # get the path to the block's parent node
                p0 = tree.find_path(s0[0], s0[1], s0[2], 3)
                if p0[3] == SectionTree.NA:
                    # The block has no sites at all
                    # Nothing to do!
                    continue

                # Now go through the block as if standard C 3d array
                for li in range(8):
                    for lj in range(8):
                        for lk in range(8):
                            # site path
                            site_path = tree.find_from(
                                s0[0] + li, s0[1] + lj, s0[2] + lk, 0, p0, 3
                            )
                            site_i = site_path[0]
                            if site_i == NA:
                                # Solid
                                bb.packer.pack_uint(0)
                            else:
                                # Fluid
                                bb.packer.pack_uint(1)
                                bb.site_count += 1
                                lnks = tree.link(site_i)
                                has_normal = False
                                if len(lnks):
                                    assert len(lnks) == 1
                                    # We have variable link data to write
                                    lnks = lnks[0]
                                    for cut_type, cut_dist, cut_bound_idx in lnks:
                                        bb.packer.pack_uint(cut_type)
                                        if cut_type == 0:
                                            continue

                                        if cut_type > 1:
                                            # IOlets need a boundary index
                                            bb.packer.pack_uint(cut_bound_idx)
                                            pass
                                        else:
                                            has_normal = True
                                        bb.packer.pack_float(cut_dist)
                                        continue

                                else:
                                    # links are all no cuts
                                    for i in range(26):
                                        bb.packer.pack_uint(0)
                                        continue
                                    pass

                                if has_normal:
                                    bb.packer.pack_uint(1)
                                    norm = tree.wall(site_i)[0]
                                    for i in range(3):
                                        bb.packer.pack_float(norm[i])
                                        pass
                                else:
                                    bb.packer.pack_uint(0)
                                    pass
                                pass

                            # Done with iter over sites


def main():
    p = argparse.ArgumentParser()
    p.add_argument("octfile", help="Input HemeLB octree file to convert to GMY")
    p.add_argument(
        "gmyfile", nargs="?", help="Output GMY file, default is basename(input) + .gmy"
    )
    args = p.parse_args()
    oct2gmy(args.octfile, args.gmyfile)


if __name__ == "__main__":
    main()
