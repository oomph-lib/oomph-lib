# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fourier_decomposed_acoustic_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fourier_decomposed_acoustic_fsi C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_unstructured_fourier_decomposed_acoustic_fsi
  SOURCES unstructured_fourier_decomposed_acoustic_fsi.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz
            oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME fourier_decomposed_acoustic_fsi
  SOURCES fourier_decomposed_acoustic_fsi.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz
            oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_fourier_decomposed_acoustic_fsi
  SOURCES unstructured_fourier_decomposed_acoustic_fsi.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz
            oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fourier_decomposed_acoustic_fsi
  TARGET_DEPENDENCIES
    adaptive_unstructured_fourier_decomposed_acoustic_fsi
    fourier_decomposed_acoustic_fsi unstructured_fourier_decomposed_acoustic_fsi
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS interaction fourier_decomposed_acoustic_fsi)

message(VERBOSE "Leaving fourier_decomposed_acoustic_fsi subdirectory")
# ------------------------------------------------------------------------------
