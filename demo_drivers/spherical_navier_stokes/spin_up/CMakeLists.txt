# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spin_up subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spin_up C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME refineable_spin_up
  SOURCES refineable_spin_up.cc
  LIBRARIES oomph::spherical_navier_stokes)

oomph_add_executable(
  NAME refineable_spin_up_cyl
  SOURCES refineable_spin_up_cyl.cc
  LIBRARIES oomph::axisym_navier_stokes)

oomph_add_executable(
  NAME spin_up
  SOURCES spin_up.cc
  LIBRARIES oomph::spherical_navier_stokes)

oomph_add_executable(
  NAME spin_up_cyl
  SOURCES spin_up_cyl.cc
  LIBRARIES oomph::axisym_navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_navier_stokes.spin_up
  TARGET_DEPENDENCIES refineable_spin_up refineable_spin_up_cyl spin_up
                      spin_up_cyl
  EXTRA_REQUIRES
  LABELS spherical_navier_stokes spin_up)

message(VERBOSE "Leaving spin_up subdirectory")
# ------------------------------------------------------------------------------
