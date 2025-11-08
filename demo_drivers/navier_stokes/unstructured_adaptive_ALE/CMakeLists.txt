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
message(VERBOSE "Entered unstructured_adaptive_ALE subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_adaptive_ALE C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME unstructured_adaptive_ALE
  SOURCES unstructured_adaptive_ALE.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::rigid_body oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_adaptive_ALE
  DEPENDS_ON unstructured_adaptive_ALE
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unstructured_adaptive_ALE subdirectory")
# ------------------------------------------------------------------------------
