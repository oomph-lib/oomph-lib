# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(prescribed_displ_lagr_mult C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult
  SOURCES prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult2
  SOURCES prescribed_displ_lagr_mult2.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult_precond
  SOURCES prescribed_displ_lagr_mult_precond.cc
  LIBRARIES oomph::multi_physics oomph::solid oomph::constitutive
            oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.prescribed_displ_lagr_mult
  TARGET_DEPENDENCIES prescribed_displ_lagr_mult prescribed_displ_lagr_mult2
                      prescribed_displ_lagr_mult_precond
  EXTRA_REQUIRES
  LABELS solid prescribed_displ_lagr_mult)

message(VERBOSE "Leaving prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
