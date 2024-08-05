# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_configuration OBJECT
        CommandLine.cc
        PathManager.cc
        SimConfig.cc
        SimBuilder.cc
        SimConfigReader.cc
        SimConfigWriter.cc
        )
