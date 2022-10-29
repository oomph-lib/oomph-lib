# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_cantilever subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(three_d_cantilever C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_cantilever
  SOURCES three_d_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME three_d_cantilever_adapt
  SOURCES three_d_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_DEFINITIONS REFINE
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.three_d_cantilever
  TARGET_DEPENDENCIES three_d_cantilever three_d_cantilever_adapt
  EXTRA_REQUIRES
  LABELS solid three_d_cantilever)

message(VERBOSE "Leaving three_d_cantilever subdirectory")
# ------------------------------------------------------------------------------
