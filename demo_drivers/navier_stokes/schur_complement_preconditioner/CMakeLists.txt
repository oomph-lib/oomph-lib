# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered schur_complement_preconditioner subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(schur_complement_preconditioner C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_fp
  SOURCES three_d_fp.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME two_d_fp
  SOURCES two_d_fp.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.schur_complement_preconditioner
  TARGET_DEPENDENCIES three_d_fp two_d_fp
  EXTRA_REQUIRES
  LABELS navier_stokes schur_complement_preconditioner)

message(VERBOSE "Leaving schur_complement_preconditioner subdirectory")
# ------------------------------------------------------------------------------
