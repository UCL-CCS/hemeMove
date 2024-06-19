#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from pathlib import Path
import os
import re
import sys
from contextlib import contextmanager


@contextmanager
def WorkingDirectory(dst_dir):
    current_dir = Path.cwd()
    os.chdir(dst_dir)
    try:
        yield dst_dir
    finally:
        os.chdir(current_dir)


# Want to match: start of line; '#include'; one or more space;
# double quote; anything (but non-greedily); double quote.
nonSystemIncludeFinder = re.compile(r'\A#include\s+"(.+?)"')
# As above, but with <...>
systemIncludeFinder = re.compile(r"\A#include\s+<(.+?)>")


def IncludedFileGenerator(filename, includeFinder):
    """Given a filename, return a generator that will yield the line
    numbers and paths included by include statements in the
    file.
    """
    with open(filename) as codefile:
        for iLine, line in enumerate(codefile):
            match = includeFinder.match(line)
            if match:
                yield (iLine + 1), match.group(1)


FileNonSystemIncludedFileGenerator = lambda filename: IncludedFileGenerator(
    filename, nonSystemIncludeFinder
)
FileSystemIncludedFileGenerator = lambda filename: IncludedFileGenerator(
    filename, systemIncludeFinder
)


cLibraryHeaders = set(
    (
        "assert.h",
        "ctype.h",
        "errno.h",
        "float.h",
        "iso646.h",
        "limits.h",
        "locale.h",
        "math.h",
        "setjmp.h",
        "signal.h",
        "stdarg.h",
        "stddef.h",
        "stdio.h",
        "stdlib.h",
        "string.h",
        "time.h",
    )
)


def CheckSystemIncludePaths(sourceFile):
    errors = False
    for lineNumber, include in FileSystemIncludedFileGenerator(sourceFile):
        if include in cLibraryHeaders:
            sys.stderr.write(
                "{file}:{line} Use of deprecated C library header <{include}>\n".format(
                    file=sourceFile, line=lineNumber, include=include
                )
            )
            errors = True
    return errors


def ignoreCopyright(f):
    line = f.readline()
    while line.strip() == "" or line[0:2] == "//":
        line = f.readline()
    return line


def GetGuardLines(filename):
    with open(filename) as f:
        lines = [ignoreCopyright(f)]
        lines.append(f.readline())

        for line in f:
            continue
        lines.append(line)
    return lines


class GuardChecker:
    def __init__(self, prefix, extensions, ignoredIncludes, skippedFiles):
        self.prefix = prefix
        self.extensions = set(extensions)
        self.ignoredIncludes = set(ignoredIncludes)
        self.skippedFiles = set(Path(f) for f in skippedFiles)

    def Walk(self, codeDir):
        for dirpath, dirnames, filenames in os.walk(codeDir):
            if "build" in dirnames:
                dirnames.remove("build")

            dirpath = Path(dirpath)
            for name in filenames:
                pth = dirpath / name
                if pth in self.skippedFiles:
                    continue
                ext = pth.suffix
                if ext in self.extensions or (
                    ext == ".in" and pth.suffixes[0] in self.extensions
                ):
                    yield pth

    def CheckIncludePaths(self, sourceFile):
        errors = False
        for lineNumber, include in FileNonSystemIncludedFileGenerator(sourceFile):
            if include in self.ignoredIncludes:
                continue

            if (
                not os.path.exists(include)
                and not os.path.exists(include + ".in")
                and not os.path.exists(os.path.splitext(include)[0] + ".in.h")
            ):
                sys.stderr.write(
                    '{file}:{line} Bad include path "{dodgy}"\n'.format(
                        file=sourceFile, line=lineNumber, dodgy=include
                    )
                )
                errors = True
        return errors

    def CheckGuardErrors(self, sourceFile):
        lines = GetGuardLines(sourceFile)

        parts = [self.prefix] + list(sourceFile.parts)
        parts[-1] = parts[-1].replace(".", "_")
        define = "_".join(parts).upper()

        error = False

        line0 = "#ifndef {define}\n".format(define=define)
        if lines[0] != line0:
            sys.stderr.write(
                "{file}:1 Bad include guard; must be {required!r} but was {actual!r}\n".format(
                    file=sourceFile, required=line0, actual=lines[0]
                )
            )
            error = True

        line1 = "#define {define}\n".format(define=define)
        if lines[1] != line1:
            sys.stderr.write(
                "{file}:2 Bad include guard; must be {required!r} but was {actual!r}\n".format(
                    file=sourceFile, required=line1, actual=lines[1]
                )
            )
            error = True

        return error

    def Check(self, directory):
        errors = False

        with WorkingDirectory(directory):
            for sourceFile in self.Walk(Path(".")):
                includeErrors = self.CheckIncludePaths(sourceFile)
                sysIncludeErrors = CheckSystemIncludePaths(sourceFile)

                ext = sourceFile.suffix
                if ext == ".h" or ext == ".hpp" or ext == ".in":
                    guardErrors = self.CheckGuardErrors(sourceFile)
                else:
                    guardErrors = False
                errors = errors or (includeErrors or sysIncludeErrors or guardErrors)

        return errors


if __name__ == "__main__":
    script = Path(sys.argv[0]).resolve()
    scriptDir = script.parent
    rootDir = scriptDir.parent
    codeDir = rootDir / "Code"

    # Main app
    codeChecker = GuardChecker(
        "HEMELB",
        (".cc", ".hpp", ".h"),
        (
            "MPWide.h",
            "tinyxml.h",
            "parmetis.h",
            "ctemplate/template.h",
            "unittests/@HEMELB_UNITTEST_INCLUDE@",
        ),
        (),
    )
    errors = codeChecker.Check(codeDir)

    # Geometry tool - common
    commonGmyChecker = GuardChecker(
        "HLBGMYTOOL_COMMON",
        (".cpp", ".hpp", ".h"),
        ("util/Vector3D.h",),
        ("PybindVTKTypeCaster.h",),
    )
    commonGmyDir = rootDir / "geometry-tool/HlbGmyTool/Model/CommonGeneration"
    errors = errors or commonGmyChecker.Check(commonGmyDir)

    # Get common includes
    with WorkingDirectory(commonGmyDir):
        common_includes = []
        for pattern in ("*.h", "*.hpp"):
            common_includes += [str(p) for p in Path(".").glob(pattern)]

    # Now old GMY
    gmyGmyChecker = GuardChecker(
        "HLBGMYTOOL_GMY",
        (".cpp", ".hpp", ".h"),
        [
            "io/formats/formats.h",
            "io/formats/geometry.h",
            "io/writers/xdr/XdrFileWriter.h",
            "io/writers/xdr/XdrMemWriter.h",
            "io/writers/xdr/XdrWriter.h",
        ]
        + common_includes,
        (
            "Point_inside_polyhedron_3.h",
            "Triangle_3_Ray_3_do_intersect.h",
        ),
    )
    gmyGmyDir = rootDir / "geometry-tool/HlbGmyTool/Model/GmyGeneration"
    errors = errors or gmyGmyChecker.Check(gmyGmyDir)

    # Last the new Octree GMY
    octGmyChecker = GuardChecker(
        "HLBGMYTOOL_OCT",
        (".cpp", ".hpp", ".h"),
        [
            "util/Vector3D.h",
            "io/formats/geometry.h",
            "Index.h",
        ]
        + common_includes,
        (),
    )
    octGmyDir = rootDir / "geometry-tool/HlbGmyTool/Model/OctGeneration"
    errors = errors or octGmyChecker.Check(octGmyDir)

    if errors:
        raise SystemExit(1)
