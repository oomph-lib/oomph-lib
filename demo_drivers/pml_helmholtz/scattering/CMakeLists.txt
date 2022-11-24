# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered scattering subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(scattering C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz_adapt
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz_scattering
  SOURCES unstructured_two_d_helmholtz_scattering.cc
  LIBRARIES oomph::pml_helmholtz)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME pml_helmholtz.scattering
  TARGET_DEPENDENCIES
    unstructured_two_d_helmholtz unstructured_two_d_helmholtz_adapt
    unstructured_two_d_helmholtz_scattering
  EXTRA_REQUIRES
  LABELS pml_helmholtz scattering)

message(VERBOSE "Leaving scattering subdirectory")
# ------------------------------------------------------------------------------
