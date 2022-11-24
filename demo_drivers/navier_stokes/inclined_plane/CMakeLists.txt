# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered inclined_plane subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(inclined_plane C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME inclined_plane
  SOURCES inclined_plane.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive)

oomph_add_executable(
  NAME inclined_plane2
  SOURCES inclined_plane.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive
  CXX_DEFINITIONS CR_ELEMENT)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.inclined_plane
  TARGET_DEPENDENCIES inclined_plane inclined_plane2
  EXTRA_REQUIRES
  LABELS navier_stokes inclined_plane)

message(VERBOSE "Leaving inclined_plane subdirectory")
# ------------------------------------------------------------------------------
