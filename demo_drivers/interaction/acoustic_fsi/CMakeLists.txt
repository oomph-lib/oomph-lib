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
message(VERBOSE "Entered acoustic_fsi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(acoustic_fsi C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME acoustic_fsi
  SOURCES acoustic_fsi.cc
  LIBRARIES oomph::helmholtz oomph::time_harmonic_linear_elasticity
            oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_acoustic_fsi
  SOURCES unstructured_acoustic_fsi.cc
  LIBRARIES oomph::helmholtz oomph::time_harmonic_linear_elasticity
            oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.acoustic_fsi
  DEPENDS_ON acoustic_fsi unstructured_acoustic_fsi
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving acoustic_fsi subdirectory")
# ------------------------------------------------------------------------------
