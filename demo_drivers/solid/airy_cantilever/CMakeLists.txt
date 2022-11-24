# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered airy_cantilever subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(airy_cantilever C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME airy_cantilever
  SOURCES airy_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive)

oomph_add_executable(
  NAME airy_cantilever2_adapt
  SOURCES airy_cantilever2.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_DEFINITIONS REFINE)

oomph_add_executable(
  NAME airy_cantilever2_noadapt
  SOURCES airy_cantilever2.cc
  LIBRARIES oomph::solid oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.airy_cantilever
  TARGET_DEPENDENCIES airy_cantilever airy_cantilever2_adapt
                      airy_cantilever2_noadapt
  EXTRA_REQUIRES
  LABELS solid airy_cantilever)

message(VERBOSE "Leaving airy_cantilever subdirectory")
# ------------------------------------------------------------------------------
