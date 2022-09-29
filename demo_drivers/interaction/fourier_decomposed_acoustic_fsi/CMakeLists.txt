# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fourier_decomposed_acoustic_fsi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fourier_decomposed_acoustic_fsi C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_unstructured_fourier_decomposed_acoustic_fsi
  SOURCES unstructured_fourier_decomposed_acoustic_fsi.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz
            oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::navier_stokes
  CXX_DEFINITIONS ADAPTIVE
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fourier_decomposed_acoustic_fsi
  SOURCES fourier_decomposed_acoustic_fsi.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz
            oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME unstructured_fourier_decomposed_acoustic_fsi
  SOURCES unstructured_fourier_decomposed_acoustic_fsi.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz
            oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fourier_decomposed_acoustic_fsi
  TARGET_DEPENDENCIES
    adaptive_unstructured_fourier_decomposed_acoustic_fsi
    fourier_decomposed_acoustic_fsi unstructured_fourier_decomposed_acoustic_fsi
  EXTRA_REQUIRES
  LABELS "interaction;fourier_decomposed_acoustic_fsi")

message(VERBOSE "Leaving fourier_decomposed_acoustic_fsi subdirectory")
# ------------------------------------------------------------------------------
