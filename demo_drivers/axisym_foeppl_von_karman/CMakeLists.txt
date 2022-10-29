# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(axisym_foeppl_von_karman C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME axisym_fvk
  SOURCES axisym_fvk.cc
  LIBRARIES oomph::axisym_foeppl_von_karman
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_foeppl_von_karman.axisym_foeppl_von_karman
  TARGET_DEPENDENCIES axisym_fvk
  EXTRA_REQUIRES
  LABELS axisym_foeppl_von_karman)

message(VERBOSE "Leaving axisym_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
