#!/usr/bin/env python
import xdrlib
import zlib
import numpy as np
from hemeTools.parsers.geometry.simple import ConfigLoader

class Decompressor(ConfigLoader):
    def __init__(self, filename, outfilename):
        ConfigLoader.__init__(self, filename)
        self.OutputFileName = outfilename
        return

    def OnEndPreamble(self):
        # Copy the preamble
        pos = self.File.tell()
        self.File.seek(0)
        
        self.OutFile = file(self.OutputFileName, 'w')
        self.OutFile.write(self.File.read(self.PreambleBytes))
        self.File.seek(pos)
        
        return
    
    def OnEndHeader(self):
        # Write dummy values
        nBlocks = self.Domain.TotalBlocks
        self.OutFile.write((4*2*nBlocks)*'\0')
        self.UncompressedBlockLengths = np.zeros(nBlocks, dtype=np.uint)
        return

    def _LoadBlock(self, domain, bIdx, bIjk):
        if domain.BlockFluidSiteCounts[bIjk] == 0:
            return
        compressed = self.File.read(self.BlockDataLength[bIjk])
        uncompressed = zlib.decompress(compressed)
        self.UncompressedBlockLengths[bIjk] = len(uncompressed)
        self.OutFile.write(uncompressed)
        return

    def OnEndBody(self):
        # Write a real header
        packer = xdrlib.Packer()
        for i in xrange(self.Domain.TotalBlocks):
            packer.pack_uint(self.Domain.BlockFluidSiteCounts[i])
            packer.pack_uint(self.UncompressedBlockLengths[i])
            
        self.OutFile.seek(self.PreambleBytes)
        self.OutFile.write(packer.get_buffer())
        self.OutFile.close()
        return
    
if __name__ == "__main__":
    import sys
    
    comp = Decompressor(sys.argv[1], sys.argv[2])
    comp.Load()
