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
message(VERBOSE "Entered extrude_quad_mesh_to_cube_mesh subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(extrude_quad_mesh_to_cube_mesh C CXX Fortran)

if(NOT oomphlib_FOUND)
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
  DEPENDS_ON extrude_triangle_generated_mesh
             extrude_with_macro_element_representation
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata triangle_meshes)

message(VERBOSE "Leaving extrude_quad_mesh_to_cube_mesh subdirectory")
# ------------------------------------------------------------------------------
