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
message(VERBOSE "Entered helmholtz_multigrid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(helmholtz_multigrid C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME structured_cubic_point_source
  SOURCES structured_cubic_point_source.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multigrid.helmholtz_multigrid
  DEPENDS_ON structured_cubic_point_source unstructured_two_d_helmholtz
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving helmholtz_multigrid subdirectory")
# ------------------------------------------------------------------------------
