# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered torus subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(torus C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME torus
  SOURCES torus.cc
  LIBRARIES oomph::navier_stokes oomph::axisym_navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.torus
  TARGET_DEPENDENCIES torus
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes torus)

message(VERBOSE "Leaving torus subdirectory")
# ------------------------------------------------------------------------------
