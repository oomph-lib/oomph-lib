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
message(VERBOSE "Entered unsteady_ring subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unsteady_ring C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME lin_unsteady_ring
  SOURCES lin_unsteady_ring.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unsteady_ring
  SOURCES unsteady_ring.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.unsteady_ring
  DEPENDS_ON lin_unsteady_ring unsteady_ring
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unsteady_ring subdirectory")
# ------------------------------------------------------------------------------
