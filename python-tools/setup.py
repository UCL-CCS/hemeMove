#!/usr/bin/env python
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from setuptools import setup, Extension
import Cython.Build
import numpy as np

ext_modules = [
    Extension(
        "hlb.utils.xdr",
        [
            "hlb/utils/xdr.pyx",
            "hlb/utils/xdr.pxd",
            "hlb/utils/cxdr.pxd",
            "hlb/utils/XdrSerialisation.c",
        ],
    ),
    Extension(
        "hlb.parsers.geometry.BaseSite",
        ["hlb/parsers/geometry/BaseSite.pyx"],
        include_dirs=[
            np.get_include(),
        ],
    ),
]

setup(
    name="hlb",
    version="0.4.0",
    description="HemeLB tools",
    author="Rupert Nash",
    author_email="r.nash@epcc.ed.ac.uk",
    packages=[
        "hlb",
        "hlb.converters",
        "hlb.parsers",
        "hlb.parsers.geometry",
        "hlb.surfacegenerator",
        "hlb.utils",
    ],
    ext_modules=ext_modules,
    cmdclass={"build_ext": Cython.Build.build_ext},
    install_requires=[
        "numpy",
        "vtk",
    ],
    tests_require=[
        "pytest",
        "pyyaml",
    ],
    entry_points={
        "console_scripts": [
            "CompressGmy = hlb.parsers.geometry.compression:compress_main",
            "DecompressGmy = hlb.parsers.geometry.compression:decompress_main",
            "DumpExtractedProperties = hlb.converters.ExtractedPropertyTextDump:main",
            "GmySelfConsistent = hlb.parsers.geometry.self_consistency:main",
            "CountSites = hlb.parsers.geometry.count_sites:main",
        ],
    },
)
