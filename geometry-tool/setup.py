# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

import os
import sys

from skbuild import setup

if sys.platform == "darwin":
    # Python thinks it's so smart and sets the
    # MACOSX_DEPLOYMENT_TARGET environment variable that messes around
    # with what features of the compiler and C++ std lib are
    # available. Set this to use your current one.
    import platform

    release, versioninfo, machine = platform.mac_ver()
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = release

setup(
    name="HemeLbSetupTool",
    version="1.2",
    author="Rupert Nash",
    author_email="r.nash@epcc.ed.ac.uk",
    packages=[
        "HemeLbSetupTool",
        "HemeLbSetupTool.Bindings",
        "HemeLbSetupTool.Util",
        "HemeLbSetupTool.Model",
        "HemeLbSetupTool.View",
        "HemeLbSetupTool.Controller",
        "HemeLbSetupTool.scripts",
    ],
    entry_points={
        "console_scripts": [
            "hemelb-gmy-cli=HemeLbSetupTool.scripts.cli:main",
            "hemelb-config2gmy=HemeLbSetupTool.scripts.config_to_geometry:main",
            "hemelb-pro2pr2=HemeLbSetupTool.scripts.pro_to_pr2:main",
        ],
        "gui_scripts": [
            "hemelb-gmy-gui=HemeLbSetupTool.scripts.gui:main[gui]",
        ],
    },
    python_requires=">=3.6",
    install_requires=[
        "pyyaml",
        # Numpy >= 1.20 requires python 3.7; VMTK conda binaries are 3.6 only
        "numpy < 1.20; python_version < '3.7'",
        "numpy; python_version >= '3.7'",
        "vtk",
        "vmtk",
    ],
    extras_require={
        "gui": ["wxPython"],
    },
)
