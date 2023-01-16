# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered harmonic subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(harmonic C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME harmonic
  SOURCES harmonic.cc
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME complex_harmonic
  SOURCES complex_harmonic.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME eigenproblems.harmonic
  TARGET_DEPENDENCIES harmonic complex_harmonic
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS eigenproblems harmonic)

message(VERBOSE "Leaving harmonic subdirectory")
# ------------------------------------------------------------------------------
