# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_free_surface subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_free_surface C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME 3d_free_surface
  SOURCES 3d_free_surface.cc
  LIBRARIES oomph::navier_stokes oomph::solid oomph::constitutive oomph::meshes
            oomph::generic oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.three_d_free_surface
  TARGET_DEPENDENCIES 3d_free_surface
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS navier_stokes three_d_free_surface)

message(VERBOSE "Leaving three_d_free_surface subdirectory")
# ------------------------------------------------------------------------------
