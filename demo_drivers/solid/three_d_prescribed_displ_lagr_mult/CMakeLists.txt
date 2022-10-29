# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(three_d_prescribed_displ_lagr_mult C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_prescribed_displ_lagr_mult
  SOURCES three_d_prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME three_d_prescribed_displ_lagr_mult_non_adapt
  SOURCES three_d_prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_DEFINITIONS NOREFINE
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.three_d_prescribed_displ_lagr_mult
  TARGET_DEPENDENCIES three_d_prescribed_displ_lagr_mult
                      three_d_prescribed_displ_lagr_mult_non_adapt
  EXTRA_REQUIRES
  LABELS solid three_d_prescribed_displ_lagr_mult)

message(VERBOSE "Leaving three_d_prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
