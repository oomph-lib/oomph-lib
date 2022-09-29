# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_linear_wave subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_linear_wave C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_linear_wave
  SOURCES two_d_linear_wave.cc
  LIBRARIES oomph::linear_wave
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_linear_wave_flux
  SOURCES two_d_linear_wave_flux.cc
  LIBRARIES oomph::linear_wave
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linear_wave.two_d_linear_wave
  TARGET_DEPENDENCIES two_d_linear_wave two_d_linear_wave_flux
  EXTRA_REQUIRES
  LABELS "linear_wave;two_d_linear_wave")

message(VERBOSE "Leaving two_d_linear_wave subdirectory")
# ------------------------------------------------------------------------------
