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
message(VERBOSE "Entered unstructured_three_d_solid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_three_d_solid C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME unstructured_three_d_solid
  SOURCES unstructured_three_d_solid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.unstructured_three_d_solid
  DEPENDS_ON unstructured_three_d_solid
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             fsi_bifurcation_solid.1.ele
             fsi_bifurcation_solid.1.face
             fsi_bifurcation_solid.1.node
             fsi_bifurcation_solid.poly
             three_d_unstructured_solid.mcr)

message(VERBOSE "Leaving unstructured_three_d_solid subdirectory")
# ------------------------------------------------------------------------------
