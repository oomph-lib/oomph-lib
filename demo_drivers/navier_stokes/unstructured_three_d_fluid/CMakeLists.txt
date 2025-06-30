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
message(VERBOSE "Entered unstructured_three_d_fluid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_three_d_fluid C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME unstructured_three_d_fluid
  SOURCES unstructured_three_d_fluid.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_three_d_fluid
  DEPENDS_ON unstructured_three_d_fluid
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             fsi_bifurcation_fluid.1.ele
             fsi_bifurcation_fluid.1.face
             fsi_bifurcation_fluid.1.node
             fsi_bifurcation_fluid.poly
             fsi_bifurcation_fluid_with_extended_tubes.1.ele
             fsi_bifurcation_fluid_with_extended_tubes.1.face
             fsi_bifurcation_fluid_with_extended_tubes.1.node
             fsi_bifurcation_fluid_with_extended_tubes.poly
             three_d_unstructured_fluid.mcr)

message(VERBOSE "Leaving unstructured_three_d_fluid subdirectory")
# ------------------------------------------------------------------------------
