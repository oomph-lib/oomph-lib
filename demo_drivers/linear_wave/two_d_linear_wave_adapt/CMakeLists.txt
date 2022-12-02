# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_linear_wave_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_linear_wave_adapt C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_linear_wave_adapt
  SOURCES two_d_linear_wave_adapt.cc
  LIBRARIES oomph::linear_wave oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linear_wave.two_d_linear_wave_adapt
  TARGET_DEPENDENCIES two_d_linear_wave_adapt
  EXTRA_REQUIRES
  LABELS linear_wave two_d_linear_wave_adapt adapt)

message(VERBOSE "Leaving two_d_linear_wave_adapt subdirectory")
# ------------------------------------------------------------------------------
