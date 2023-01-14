# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered extrude_quad_mesh_to_cube_mesh subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(extrude_quad_mesh_to_cube_mesh C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME extrude_triangle_generated_mesh
  SOURCES extrude_triangle_generated_mesh.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME extrude_with_macro_element_representation
  SOURCES extrude_with_macro_element_representation.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME extrude_quad_mesh_to_cube_mesh.extrude_quad_mesh_to_cube_mesh
  TARGET_DEPENDENCIES extrude_triangle_generated_mesh
                      extrude_with_macro_element_representation
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES triangle_meshes
  LABELS extrude_quad_mesh_to_cube_mesh)

message(VERBOSE "Leaving extrude_quad_mesh_to_cube_mesh subdirectory")
# ------------------------------------------------------------------------------
