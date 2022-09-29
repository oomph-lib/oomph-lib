# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vmtk_fluid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(vmtk_fluid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME vmtk_fluid
  SOURCES vmtk_fluid.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.vmtk_fluid
  TARGET_DEPENDENCIES vmtk_fluid
  EXTRA_REQUIRES fluid_iliac.1.ele fluid_iliac.1.face fluid_iliac.1.node
                 fluid_iliac.poly iliac.xda
  LABELS "navier_stokes;vmtk_fluid")

message(VERBOSE "Leaving vmtk_fluid subdirectory")
# ------------------------------------------------------------------------------
