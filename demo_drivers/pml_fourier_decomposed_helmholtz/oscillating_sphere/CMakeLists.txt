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
message(VERBOSE "Entered oscillating_sphere subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(oscillating_sphere C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adaptive_oscillating_sphere
  SOURCES oscillating_sphere.cc
  LIBRARIES oomph::pml_fourier_decomposed_helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME oscillating_sphere
  SOURCES oscillating_sphere.cc
  LIBRARIES oomph::pml_fourier_decomposed_helmholtz oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME pml_fourier_decomposed_helmholtz.oscillating_sphere
  DEPENDS_ON adaptive_oscillating_sphere oscillating_sphere
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving oscillating_sphere subdirectory")
# ------------------------------------------------------------------------------
