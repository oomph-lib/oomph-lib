# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered quad_from_triangle_mesh subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(quad_from_triangle_mesh C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_moving_block_navier_stokes
  SOURCES adaptive_moving_block_navier_stokes.cc
  LIBRARIES oomph::navier_stokes oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME adaptive_unstructured_scattering_quad
  SOURCES adaptive_unstructured_scattering_quad.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME adaptive_unstructured_two_d_solid
  SOURCES adaptive_unstructured_two_d_solid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_two_d_solid
  SOURCES unstructured_two_d_solid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.quad_from_triangle_mesh
  TARGET_DEPENDENCIES
    adaptive_moving_block_navier_stokes adaptive_unstructured_scattering_quad
    adaptive_unstructured_two_d_solid unstructured_two_d_solid
  EXTRA_REQUIRES
    box_hole.1.ele
    box_hole.1.node
    box_hole.1.poly
    box_hole.poly
    solid.1.ele
    solid.1.node
    solid.1.poly
    solid.poly
  LABELS meshing quad_from_triangle_mesh)

message(VERBOSE "Leaving quad_from_triangle_mesh subdirectory")
# ------------------------------------------------------------------------------
