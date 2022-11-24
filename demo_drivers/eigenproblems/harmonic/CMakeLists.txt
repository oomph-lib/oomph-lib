# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered harmonic subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(harmonic C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

if(HAVE_ARPACK_SOURCES)
  oomph_add_executable(
    NAME harmonic
    SOURCES harmonic.cc
    LIBRARIES oomph::generic)

  oomph_add_executable(
    NAME complex_harmonic
    SOURCES complex_harmonic.cc
    LIBRARIES oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS)

if(HAVE_ARPACK_SOURCES)
  list(APPEND TEST_TARGET_DEPS harmonic complex_harmonic)
endif()

oomph_add_test(
  TEST_NAME eigenproblems.harmonic
  TARGET_DEPENDENCIES ${TEST_TARGET_DEPS}
  EXTRA_REQUIRES
  LABELS eigenproblems harmonic)

message(VERBOSE "Leaving harmonic subdirectory")
# ------------------------------------------------------------------------------
