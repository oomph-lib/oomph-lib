# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered thermo subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(thermo C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME thermo
  SOURCES thermo.cc
  LIBRARIES oomph::unsteady_heat oomph::solid oomph::constitutive oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.thermo
  TARGET_DEPENDENCIES thermo
  EXTRA_REQUIRES
  LABELS multi_physics thermo)

message(VERBOSE "Leaving thermo subdirectory")
# ------------------------------------------------------------------------------
