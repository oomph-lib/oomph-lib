# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(extrude_quad_mesh_to_cube_mesh C CXX Fortran)
find_package(oomphlib REQUIRED)
include(CTest)

oomph_add_executable(
  NAME extrude_triangle_generated_mesh
  SOURCES extrude_triangle_generated_mesh.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME extrude_with_macro_element_representation
  SOURCES extrude_with_macro_element_representation.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME extrude_quad_mesh_to_cube_mesh.extrude_quad_mesh_to_cube_mesh
  TARGET_DEPENDENCIES extrude_triangle_generated_mesh
                      extrude_with_macro_element_representation
  EXTRA_REQUIRES triangle_meshes
  LABELS "extrude_quad_mesh_to_cube_mesh")

# ------------------------------------------------------------------------------
