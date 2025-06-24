# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered cylinder subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(cylinder C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adaptive_pressure_loaded_cylinder
  SOURCES pressure_loaded_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME adaptive_unstructured_cylinder
  SOURCES unstructured_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME cylinder
  SOURCES cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME pressure_loaded_cylinder
  SOURCES pressure_loaded_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_cylinder
  SOURCES unstructured_cylinder.cc
  LIBRARIES oomph::time_harmonic_fourier_decomposed_linear_elasticity
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME time_harmonic_fourier_decomposed_linear_elasticity.cylinder
  DEPENDS_ON adaptive_pressure_loaded_cylinder adaptive_unstructured_cylinder
             cylinder pressure_loaded_cylinder unstructured_cylinder
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving cylinder subdirectory")
# ------------------------------------------------------------------------------
