# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

find_hemelb_dependency(Catch2 REQUIRED)

# Sources files in ARGN
function(add_test_executable target)
  add_executable(${target} ${ARGN})
  target_compile_definitions(${target} PRIVATE HEMELB_DOING_UNITTESTS)
  # ReporterTests directly use ctemplate
  # FlowExtensionTests directly uses tinyxml
  target_link_libraries(${target} PRIVATE
    Catch2::Catch2
    ${heme_libraries}
    ${MPI_LIBRARIES}
    test-helpers
    )
  if (HEMELB_BUILD_RBC)
    # redblood tests need VTK
    target_link_libraries(${target} PRIVATE VTK::CommonCore)
  endif()

  codesign(${target})
  install(TARGETS ${target} RUNTIME DESTINATION bin)
endfunction()

# Source files in ARGN
function(add_test_lib name)
  add_library(${name} OBJECT ${ARGN})
  target_link_libraries(${name} PRIVATE Catch2::Catch2)
  target_compile_definitions(${name} PUBLIC HEMELB_DOING_UNITTESTS)
endfunction()

add_to_resources(resources/four_cube.gmy resources/four_cube.xml resources/four_cube_multiscale.xml
  resources/config.xml resources/config0_2_0.xml
  resources/config_file_inlet.xml resources/iolet.txt
  resources/config-velocity-iolet.xml resources/config_new_velocity_inlets.xml
  resources/velocity_inlet.txt.weights.txt
  resources/xmltest.xml resources/config_file_velocity_inlet.xml resources/velocity_inlet.txt
  resources/large_cylinder.gmy resources/large_cylinder.xml
)


add_test_executable(hemelb-tests main.cc SimulationControllerTests.cc)

add_subdirectory(helpers)

macro(add_to_tests dirname)
  add_subdirectory(${dirname})
  target_link_libraries(hemelb-tests PRIVATE "test_${dirname}")
endmacro()

add_to_tests(configuration)
add_to_tests(extraction)
add_to_tests(geometry)
add_to_tests(io)
add_to_tests(lb)
add_to_tests(net)
add_to_tests(reporting)
add_to_tests(util)

if (HEMELB_BUILD_MULTISCALE)
  add_to_resources(multiscale/mpwide/MPWSettings.cfg)
  add_to_tests(multiscale)
endif()

if (HEMELB_BUILD_RBC)
  add_to_resources(
    resources/red_blood_cell.txt resources/red_blood_cube.txt
    resources/large_cylinder_rbc.xml
    resources/rbc_ico_1280.msh resources/rbc_ico_720.msh
    resources/fedosov1c.xml resources/fedosov1c.gmy
    resources/sad.msh
    resources/cyl_l100_r5.xml
    resources/cyl_l100_r5.gmy
    resources/rbc_ico_2880.msh
    resources/rbc_ico_720_correct.msh
    resources/rbc_ico_720.msh
    resources/rbc_ico_720.vtp
    resources/992Particles_rank3_26_t992.msh
    resources/992Particles_rank3_26_t992.vtp
    resources/empty_for_relative_paths.xml
    )

  add_to_tests(redblood)
  add_test_executable(mpi_redblood_tests main.cc)
  target_link_libraries(mpi_redblood_tests PRIVATE test_mpi_redblood)
  add_executable(hemelb-tetgen redblood/tetrahedron.cc)
  target_link_libraries(hemelb-tetgen ${heme_libraries} ${MPI_LIBRARIES})
endif()
