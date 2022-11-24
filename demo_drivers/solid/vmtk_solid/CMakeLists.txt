# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vmtk_solid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vmtk_solid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME vmtk_solid
  SOURCES vmtk_solid.cc
  LIBRARIES oomph::solid oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.vmtk_solid
  TARGET_DEPENDENCIES vmtk_solid
  EXTRA_REQUIRES solid_iliac.1.ele solid_iliac.1.face solid_iliac.1.node
                 solid_iliac.poly
  LABELS solid vmtk_solid)

message(VERBOSE "Leaving vmtk_solid subdirectory")
# ------------------------------------------------------------------------------
