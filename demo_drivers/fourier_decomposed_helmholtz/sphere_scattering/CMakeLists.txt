# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered sphere_scattering subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(sphere_scattering C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME sphere_scattering
  SOURCES sphere_scattering.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_adaptive_sphere_scattering
  SOURCES unstructured_sphere_scattering.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_sphere_scattering
  SOURCES unstructured_sphere_scattering.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME fourier_decomposed_helmholtz.sphere_scattering
  TARGET_DEPENDENCIES sphere_scattering unstructured_adaptive_sphere_scattering
                      unstructured_sphere_scattering
  EXTRA_REQUIRES
  LABELS fourier_decomposed_helmholtz sphere_scattering)

message(VERBOSE "Leaving sphere_scattering subdirectory")
# ------------------------------------------------------------------------------
