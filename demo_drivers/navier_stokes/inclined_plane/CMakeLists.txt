# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered inclined_plane subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(inclined_plane C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME inclined_plane
  SOURCES inclined_plane.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME inclined_plane2
  SOURCES inclined_plane.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic
  CXX_DEFINITIONS CR_ELEMENT)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.inclined_plane
  DEPENDS_ON inclined_plane inclined_plane2
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving inclined_plane subdirectory")
# ------------------------------------------------------------------------------
