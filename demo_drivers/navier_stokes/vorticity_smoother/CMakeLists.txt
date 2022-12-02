# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vorticity_smoother subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vorticity_smoother C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vorticity_smoother_validation
  SOURCES vorticity_smoother_validation.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.vorticity_smoother
  TARGET_DEPENDENCIES vorticity_smoother_validation
  EXTRA_REQUIRES
  LABELS navier_stokes vorticity_smoother)

message(VERBOSE "Leaving vorticity_smoother subdirectory")
# ------------------------------------------------------------------------------
