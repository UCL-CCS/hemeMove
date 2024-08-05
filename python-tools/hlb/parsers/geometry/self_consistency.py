#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import collections
import os.path
import queue
import threading

import numpy as np

from . import MooreNeighbourhoodDirections
from .generic import (
    OutOfDomainSite,
    AllSolidBlock,
    Site,
    Block,
)
from .multiprocess import AsyncBlockProcessingLoader

REALISTIC_IOLET_COUNT = 100


class EmptyErrorCollectionFormatError(Exception):
    pass


class Error(object):
    """Base class to hold info about the error in a config file."""

    def __init__(self, message):
        self.message = message
        return

    _template = "{message}"

    def Format(self):
        return self._template.format(**self.__dict__)

    def Indent(self, text):
        return "    " + text

    pass


class DomainError(Error):
    """Config file errors to do with domain global data."""

    pass


class ErrorCollection(Error):
    """A container for Errors."""

    def __init__(self, item):
        self.item = item
        self.itemErrors = []
        self.subItemErrors = collections.OrderedDict()
        return

    @property
    def TotalErrors(self):
        return len(self.itemErrors) + len(self.subItemErrors)

    def Format(self):
        nErr = self.TotalErrors
        if nErr == 0:
            return []
        elif nErr == 1:
            self.errString = "Error"
        else:
            self.errString = "Errors"
            pass

        errLines = ["{errString} in {item}:".format(**self.__dict__)]
        for iErr in self.itemErrors:
            errLines.append(self.Indent(iErr.Format()))
            continue

        for siErrColl in self.subItemErrors:
            try:
                subLines = siErrColl.Format()
                errLines.extend(self.Indent(sl) for sl in subLines)
            except EmptyErrorCollectionFormatError:
                pass
            continue

        return errLines

    pass


class BlockErrorCollection(ErrorCollection):
    """Container for errors to do with a particular block."""

    def AddBlockError(self, message):
        self.itemErrors.append(BlockError(self.item, message))
        return

    def AddSiteError(self, site, message):
        try:
            siteErrorCollection = self.subItemErrors[site]
        except KeyError:
            self.subItemErrors[site] = SiteErrorCollection(site)
            siteErrorCollection = self.subItemErrors[site]
            pass

        siteErrorCollection.Add(message)
        return

    pass


class SiteErrorCollection(ErrorCollection):
    """Container for errors to do with a particular site."""

    def Add(self, message):
        self.itemErrors.append(SiteError(self.item, message))
        return

    pass


class BlockError(Error):
    """An error attached a block's data."""

    _template = "{block}: {message}"

    def __init__(self, block, message):
        Error.__init__(self, message)
        self.block = block
        return

    pass


class SiteError(Error):
    """An error in a site's data."""

    _template = "{message}"

    def __init__(self, site, message):
        Error.__init__(self, message)
        self.site = site
        self.block = site.GetBlock()
        return

    pass


class CheckingLoader(AsyncBlockProcessingLoader):
    BLOCK_REPORT_PERIOD = 100

    def __init__(self, *args, **kwargs):
        self.Verbose = kwargs.pop("verbose", False)
        AsyncBlockProcessingLoader.__init__(self, *args, **kwargs)
        return

    def Info(self, message):
        if self.Verbose:
            print(message)
        return

    def Load(self):
        self.StartReportingThread()
        # Use a try finally to ensure the error thread gets cleaned
        # up, even if we get an exception.
        try:
            ans = AsyncBlockProcessingLoader.Load(self)
        finally:
            self.EndReportingThread()

        return ans

    def OnBlockProcessed(self, bIdx):
        bIjk = self.Domain.BlockIndexer.NdToOne(bIdx)
        if bIjk % self.BLOCK_REPORT_PERIOD == 0:
            self.Info("Checked block {} of {}".format(bIjk, self.Domain.TotalBlocks))
            pass
        return AsyncBlockProcessingLoader.OnBlockProcessed(self, bIdx)

    def OnEndPreamble(self):
        """Loaded the very basic domain information. Give a status
        report.
        """
        fields = ["Version", "BlockSize", "BlockCounts", "VoxelSize", "Origin"]
        template = "\n".join("%s: {0.%s}" % (f, f) for f in fields)
        self.Info(template.format(self.Domain))
        return

    def OnEndHeader(self):
        """All header data now loaded. Give a status report and check
        consistency of header data with itself and the actual file
        size.
        """
        fluidSiteCount = np.sum(self.Domain.BlockFluidSiteCounts)
        self.Info("NumberOfFluidSites: {}".format(fluidSiteCount))

        # For consistency, if BlockDataLength[i] == 0 then
        # BlockFluidSiteCounts[i] must also be zero, and vice versa
        for bIjk, bIdx in self.Domain.BlockIndexer.IterBoth():
            if (
                self.BlockDataLength[bIjk] == 0
                and self.Domain.BlockFluidSiteCounts[bIjk] != 0
            ):
                self.PrintError(
                    BlockError(
                        self.Domain.GetBlock(bIdx),
                        "Header states no data but specifies some " "fluid sites",
                    ).Format()
                )
                pass

            if (
                self.Domain.BlockFluidSiteCounts[bIjk] == 0
                and self.BlockDataLength[bIjk] != 0
            ):
                self.PrintError(
                    BlockError(
                        self.Domain.GetBlock(bIdx),
                        "Header states no fluid sites but specifies " "some data",
                    ).Format()
                )
                pass
            continue

        # The length of the file must be equal to the value we
        # calculate from the headers.
        claimedFileSize = (
            self.PreambleBytes + self.HeaderBytes + np.sum(self.BlockDataLength)
        )
        if claimedFileSize != os.path.getsize(self.GmyFileName):
            self.PrintError(
                DomainError("File length does not match file metadata").Format()
            )
            pass

        self.Checker = BlockChecker()
        self.SetBlockProcessor(self.Checker)
        return

    def HandleBlockProcessingResult(self, result):
        # Result is a list of errors; only bother sending it if there
        # are items
        if len(result):
            self.ReportQueue.put(result)
        return

    def PrintError(self, err):
        self.ReportQueue.put([err])
        return

    def StartReportingThread(self):
        """Start a thread to deal with writing our reports."""
        self.ReportQueue = queue.Queue()
        self.ReportingThread = threading.Thread(
            target=self.ReportingThreadMain, name="ReportingThread"
        )
        self.ReportingThread.start()
        return

    def ReportingThreadMain(self):
        """This is the entry point for the reporting thread. Make sure
        we're thread-safe here!

        We will just take items off the queue and print them unless
        they're the sentinal value (None), in which case we return,
        killing the reporing thread.
        """
        hadErrors = False
        while True:
            # ReportQueue's a Queue.Queue so this is OK
            # (It will block until we get something)
            errList = self.ReportQueue.get()
            if errList is None:
                # Got the sentinal, so we need to exit
                # First, set a flag to indicate the presence of errors.

                # This assignment is OK since the main thread waits on
                # us after putting the sentinal into the queue (see
                # EndReportingThread).
                self.HadErrors = hadErrors
                break
            hadErrors = True
            for e in errList:
                print(e)

    def EndReportingThread(self):
        """Kill the reporting thread by sending the sentinal value and
        clean up.
        """
        self.ReportQueue.put(None)
        self.ReportingThread.join()
        return

    pass


class BlockChecker(object):
    """Callable object that is run in a separate process to check a
    single block for consistency.

    """

    def __call__(self, domain, bInd):
        """Now we iterate over every site on the block, and check
        constraints as specified in the Check* methods.

        A list of error STRINGS is returned, i.e. the result of
        calling Format() on the accumulated Error object.
        """

        # So this may LOOK dodgy, using an instance attribute to store
        # the errors when we've only instantiated one of these and are
        # checking a load of blocks in parallel. HOWEVER, we are
        # pickling the unexecuted instance created above out to each
        # process to be run per-block. Once results have been sent
        # back, the executed object is destroyed.
        block = domain.GetBlock(bInd)
        self.block = block
        blockErrors = BlockErrorCollection(block)
        self.numFluid = 0

        dom = block.GetDomain()

        if isinstance(block, AllSolidBlock):
            self.CheckFluidSiteCount(blockErrors.AddBlockError)
            return blockErrors.Format()

        # gs = global site
        # ls = local site
        for lsInd in np.ndindex((dom.BlockSize, dom.BlockSize, dom.BlockSize)):
            lsInd = np.array(lsInd)

            site = block.GetLocalSite(lsInd)

            # Quick helper function to add an error to this site's
            # error list
            addSiteError = lambda msg: blockErrors.AddSiteError(site, msg)

            self.CheckBasicFluidProperty(site, addSiteError)

            if not site.IsFluid:
                # Solid sites we don't check in detail since they have
                # no more data
                continue

            self.numFluid += 1

            self.CheckFluidSiteLinks(site, addSiteError)

            self.CheckFluidSiteWallNormal(site, addSiteError)

            continue

        self.CheckFluidSiteCount(blockErrors.AddBlockError)
        return blockErrors.Format()

    def CheckFluidSiteCount(self, addBlockError):
        """The number of non-solid sites on each block must equal that
        declared in the header.
        """
        if self.numFluid != self.block.nFluidSites:
            addBlockError("File fluid site counts not equal to number loaded")
            pass
        return

    def CheckBasicFluidProperty(self, site, addSiteError):
        """Some basic sanity checking on the value of the 'isFluid' uint."""
        if not (site.IsFluid == Site.SOLID or site.IsFluid == Site.FLUID):
            addSiteError(
                "Site doesn't appear to have an appropriate value for 'IsFluid': %d".format(
                    site.IsFluid
                )
            )
        return

    def CheckFluidSiteLinks(self, site, addSiteError):
        """Loop over the Moore neighbourhood, checking link properties
        and whether the type is consistent with the neighbours' types.
        More concretely:

        - links must have a valid intersection type;

        - links to solid sites must have the intersection type set to
          not NO_INTERSECTION

        - links to fluid sites must have intersection type
          NO_INTERSECTION;

        - links with an intersection must have a cut distance in (0,1);

        - links with type INLET/OUTLET must have a valid IOlet index; and

        - links with non-IOlet type must not have a valid IOlet index.
        """
        for iNeigh, delta in enumerate(MooreNeighbourhoodDirections):
            neigh = site.GetBlock().GetSite(site.Index + delta)
            linkType = site.IntersectionType[iNeigh]
            # Check link type is one we know about
            if linkType not in Site.INTERSECTION_TYPES:
                addSiteError(
                    "Site had an invalid intersection type: {0}".format(linkType)
                )

            # Check site type consistency with neighbour
            if isinstance(neigh, OutOfDomainSite) or neigh.IsSolid:
                # neigh is solid
                if linkType == Site.NO_INTERSECTION:
                    addSiteError(
                        "Site has no intersection on link to solid neighbour {0}".format(
                            neigh
                        )
                    )
                pass
            else:
                # neigh is fluid
                if linkType != Site.NO_INTERSECTION:
                    addSiteError(
                        "Site has intersection (type {0}) on link to fluid neighbour {1}".format(
                            linkType, neigh
                        )
                    )
                    pass
                pass

            # If intersection, check distance in (0,1)
            if linkType != Site.NO_INTERSECTION:
                distance = site.IntersectionDistance[iNeigh]
                if distance < 0.0 or distance >= 1.0:
                    addSiteError(
                        "Site had an intersection distance (%f) outside the allowed range (0,1)"
                        % distance
                    )
                    pass
                pass

            # If IOlet, check the index
            ioletIndex = site.IOletIndex[iNeigh]
            if linkType in [Site.INLET_INTERSECTION, Site.OUTLET_INTERSECTION]:
                if ioletIndex < 0 or ioletIndex > REALISTIC_IOLET_COUNT:
                    addSiteError(
                        "IOlet index for an iolet link (%d) was outside "
                        "expected range.".format(ioletIndex)
                    )
                    pass
            elif ioletIndex != Site.NO_IOLET:
                addSiteError(
                    "Site had no link crossing an iolet but DID have "
                    "an iolet index set (%d).".format(ioletIndex)
                )
                pass
            continue
        return

    def CheckFluidSiteWallNormal(self, site, addSiteError):
        """Check that a wall normal is available for sites with a link
        intersecting a wall. Check that the provided normal is a unit
        vector.
        """
        if site.WallNormalAvailable:
            if not np.any(site.IntersectionType == Site.WALL_INTERSECTION):
                addSiteError(
                    "Site claims to have a wall normal available but no link intersects a wall"
                )
            if np.fabs(np.linalg.norm(site.WallNormal) - 1.0) > 1e-6:
                addSiteError("Site has a non unitarian wall normal")
        else:
            if np.any(site.IntersectionType == Site.WALL_INTERSECTION):
                addSiteError(
                    "Site has a link intersecting a wall but no wall normal available"
                )
        return

    pass


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze a HemeLB input file for self-consistency.",
        epilog="The return value indicates the presence of detected errors.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=False,
        help="enable quiet mode, which supresses non-error output",
    )
    parser.add_argument(
        "--debug", action="store_true", default=False, help="enable debugging output"
    )
    parser.add_argument("input", nargs=1, help="the input file to check")
    args = parser.parse_args()

    if args.debug:
        from .multiprocess import GetLogger

        GetLogger().setLevel("DEBUG")
        pass

    ldr = CheckingLoader(args.input[0], verbose=(not args.quiet))
    ldr.Load()

    if ldr.HadErrors:
        raise SystemExit(1)
