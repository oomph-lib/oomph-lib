# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered elastic_bretherton subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(elastic_bretherton C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_bretherton
  SOURCES elastic_bretherton.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::fluid_interface
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.elastic_bretherton
  TARGET_DEPENDENCIES elastic_bretherton
  EXTRA_REQUIRES
  LABELS "interaction;elastic_bretherton")

message(VERBOSE "Leaving elastic_bretherton subdirectory")
# ------------------------------------------------------------------------------
