# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pseudo_solid_collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(pseudo_solid_collapsible_channel C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME pseudo_solid_collapsible_channel
  SOURCES pseudo_solid_collapsible_channel.cc
  LIBRARIES oomph::navier_stokes oomph::solid oomph::constitutive
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.pseudo_solid_collapsible_channel
  TARGET_DEPENDENCIES pseudo_solid_collapsible_channel
  EXTRA_REQUIRES
  LABELS "solid;pseudo_solid_collapsible_channel")

message(VERBOSE "Leaving pseudo_solid_collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
