# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

function(get_vtkIdType outvar)
  # Figure out what vtkIdType is
  list(APPEND CMAKE_REQUIRED_LIBRARIES VTK::CommonCore)

  # These are the possibilities
  set(possible_idtypes "long long" long int)
  foreach(potential_idtype ${possible_idtypes})
    message(STATUS "idtype=${potential_idtype}")
    # Spaces in variable names cause problems
    string(REPLACE " " "_" varname "${potential_idtype}")
    string(TOUPPER ${varname} varname)

    check_cxx_source_compiles("
#include <type_traits>
#include <vtkType.h>
static_assert(std::is_same<vtkIdType, ${potential_idtype}>::value, \"\");
int main(){}
"
      HAVE_VTKIDTYPE_${varname})

    if (HAVE_VTKIDTYPE_${varname})
      set(${outvar} "${potential_idtype}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  message(FATAL_ERROR "Could not figure out vtkIdType")
endfunction()

if (HEMELB_BUILD_RBC)
  # Build full RBC code
  get_vtkIdType(HEMELB_RBC_IDTYPE)
  configure_file (
  "${PROJECT_SOURCE_DIR}/redblood/MeshIdType.h.in"
  "${PROJECT_BINARY_DIR}/redblood/MeshIdType.h"
  )

  add_library(hemelb_redblood OBJECT
    CellControllerBuilder.cc
    Mesh.cc MeshIO.cc
    CellIO.cc
    CellBase.cc Cell.cc CellEnergy.cc Facet.cc
    Interpolation.cc
    CellCell.cc FlowExtension.cc FaderCell.cc RBCInserter.cc
    VertexBag.cc Borders.cc
    WallCellPairIterator.cc
    VTKError.cc
    buffer/Buffer.cc buffer/Columns.cc
    parallel/NodeCharacterizer.cc parallel/CellParallelization.cc parallel/SpreadForces.cc
    parallel/IntegrateVelocities.cc parallel/GraphBasedCommunication.cc
    )

  target_link_libraries(hemelb_redblood PRIVATE
    hdf5::hdf5_hl
    VTK::CommonCore VTK::CommonDataModel VTK::IOXML VTK::FiltersCore)
else()
  # Build minimal interface
  add_library(hemelb_redblood OBJECT
    FlowExtension.cc
    )
endif()
