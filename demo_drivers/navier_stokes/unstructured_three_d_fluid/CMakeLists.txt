# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_three_d_fluid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
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
             three_d_unstructured_fluid.mcr
  LABELS navier_stokes unstructured_three_d_fluid)

message(VERBOSE "Leaving unstructured_three_d_fluid subdirectory")
# ------------------------------------------------------------------------------
