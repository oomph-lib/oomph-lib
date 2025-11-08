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
message(VERBOSE "Entered unstructured_adaptive_fs subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_adaptive_fs C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adaptive_bubble_in_channel
  SOURCES adaptive_bubble_in_channel.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::fluid_interface oomph::meshes oomph::generic)

oomph_add_executable(
  NAME adaptive_drop_in_channel
  SOURCES adaptive_drop_in_channel.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::fluid_interface oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_adaptive_fs
  DEPENDS_ON adaptive_bubble_in_channel adaptive_drop_in_channel
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unstructured_adaptive_fs subdirectory")
# ------------------------------------------------------------------------------
