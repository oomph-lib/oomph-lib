# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered orr_sommerfeld subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(orr_sommerfeld C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

if(HAVE_ARPACK_SOURCES)
  oomph_add_executable(
    NAME orr_sommerfeld
    SOURCES orr_sommerfeld.cc
    LIBRARIES oomph::generic
    CXX_STANDARD 14)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS)

if(HAVE_ARPACK_SOURCES)
  list(APPEND TEST_TARGET_DEPS orr_sommerfeld)
endif()

oomph_add_test(
  TEST_NAME eigenproblems.orr_sommerfeld
  TARGET_DEPENDENCIES ${TEST_TARGET_DEPS}
  EXTRA_REQUIRES
  LABELS "eigenproblems;orr_sommerfeld")

message(VERBOSE "Leaving orr_sommerfeld subdirectory")
# ------------------------------------------------------------------------------
