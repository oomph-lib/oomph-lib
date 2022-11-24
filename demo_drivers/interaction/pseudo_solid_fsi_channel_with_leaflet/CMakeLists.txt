# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pseudo_solid_fsi_channel_with_leaflet subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(pseudo_solid_fsi_channel_with_leaflet C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_channel_with_leaflet_precond
  SOURCES fsi_channel_with_leaflet_precond.cc
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes
            oomph::constitutive oomph::solid)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.pseudo_solid_fsi_channel_with_leaflet
  TARGET_DEPENDENCIES fsi_channel_with_leaflet_precond
  EXTRA_REQUIRES
  LABELS interaction pseudo_solid_fsi_channel_with_leaflet)

message(VERBOSE "Leaving pseudo_solid_fsi_channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
