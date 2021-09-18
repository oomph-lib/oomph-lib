# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_vmtk_fsi subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(unsteady_vmtk_fsi C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME unsteady_vmtk_fsi
  SOURCES unsteady_vmtk_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::poisson
            oomph::linear_elasticity
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unsteady_vmtk_fsi
  TARGET_DEPENDENCIES unstructured_adaptive_2d_fsi
  EXTRA_REQUIRES
    fluid.1.ele
    fluid.1.face
    fluid.1.node
    fluid_quadratic_nodes.dat
    solid.1.ele
    solid.1.face
    solid.1.node
    solid_quadratic_nodes.dat
    boundary_enumeration.dat
  LABELS "interaction;unsteady_vmtk_fsi;adapt")

message(VERBOSE "Leaving unsteady_vmtk_fsi subdirectory")
# ------------------------------------------------------------------------------
