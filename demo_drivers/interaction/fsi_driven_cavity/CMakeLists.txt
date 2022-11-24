# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_driven_cavity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_driven_cavity_driver
  SOURCES fsi_driven_cavity_driver.cc fsi_driven_cavity_problem.h
  LIBRARIES oomph::beam oomph::navier_stokes)

oomph_add_executable(
  NAME fsi_driven_cavity_segregated_driver
  SOURCES fsi_driven_cavity_segregated_driver.cc fsi_driven_cavity_problem.h
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_driven_cavity
  TARGET_DEPENDENCIES fsi_driven_cavity_driver
                      fsi_driven_cavity_segregated_driver
  EXTRA_REQUIRES
  LABELS interaction fsi_driven_cavity
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving fsi_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
