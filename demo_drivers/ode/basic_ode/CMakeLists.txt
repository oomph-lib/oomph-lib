# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered basic_ode subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(basic_ode C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME basic_ode
  SOURCES basic_ode.cc
  LIBRARIES oomph::ode oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME ode.basic_ode
  TARGET_DEPENDENCIES basic_ode
  EXTRA_REQUIRES zeros
  LABELS ode basic_ode
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving basic_ode subdirectory")
# ------------------------------------------------------------------------------
