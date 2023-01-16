# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_interface subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(adaptive_interface C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_interface
  SOURCES adaptive_interface.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.adaptive_interface
  TARGET_DEPENDENCIES adaptive_interface
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS navier_stokes adaptive_interface adapt)

message(VERBOSE "Leaving adaptive_interface subdirectory")
# ------------------------------------------------------------------------------
