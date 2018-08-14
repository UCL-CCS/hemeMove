hemelb_dependency(MPWide find)
if (MPWide_FOUND)
  message("MPWide already installed, no need to (download and) build")
  add_custom_target(MPWide)
else()
  message("MPWide not installed, will build from source")
  find_file(MPWIDE_TARBALL MPWide-1.2b.tar.gz 
  DOC "Path to download MPWide (can be url http://)"
  PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
  )
if(NOT MPWIDE_TARBALL)
  message("No MPWide source found, will download.")
  set(MPWIDE_TARBALL http://castle.strw.leidenuniv.nl/~derek/MPWide-1.2b.tar.gz
    CACHE STRING "Path to download MPWide (can be local file://)" FORCE)
endif()
set(PATCH_COMMAND_MPWIDE patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/mpwide_include.diff )
ExternalProject_Add(
  MPWide
  URL ${MPWIDE_TARBALL}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND echo "No configure step for MPWIDE"
  BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS} 
  PATCH_COMMAND ${PATCH_COMMAND_MPWIDE}
  INSTALL_COMMAND cp <SOURCE_DIR>/libMPW.a ${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib && cp <SOURCE_DIR>/MPWide.h ${HEMELB_DEPENDENCIES_INSTALL_PATH}/include
  )
endif()
