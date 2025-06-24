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
message(VERBOSE "Entered unstructured_two_d_curved subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_two_d_curved C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adaptive_unstructured_two_d_curved
  SOURCES unstructured_two_d_curved.cc
  LIBRARIES oomph::axisym_poroelasticity oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_two_d_curved
  SOURCES unstructured_two_d_curved.cc
  LIBRARIES oomph::axisym_poroelasticity oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_poroelasticity.unstructured_two_d_curved
  DEPENDS_ON adaptive_unstructured_two_d_curved unstructured_two_d_curved
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unstructured_two_d_curved subdirectory")
# ------------------------------------------------------------------------------
