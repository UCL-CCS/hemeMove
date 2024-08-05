# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_extraction OBJECT
  GeometrySelector.cc
  StraightLineGeometrySelector.cc LocalPropertyOutput.cc
  PlaneGeometrySelector.cc PropertyActor.cc
  PropertyWriter.cc WholeGeometrySelector.cc LbDataSourceIterator.cc
  GeometrySurfaceSelector.cc SurfacePointSelector.cc LocalDistributionInput.cc)
