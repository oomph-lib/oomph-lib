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
message(VERBOSE "Entered mesh_from_inline_triangle subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(mesh_from_inline_triangle C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME mesh_from_inline_triangle
  SOURCES mesh_from_inline_triangle.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_inline_triangle_no_adapt
  SOURCES mesh_from_inline_triangle_no_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_inline_triangle_polygon
  SOURCES mesh_from_inline_triangle.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic
  CXX_DEFINITIONS OUTER_POLYGON)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_inline_triangle
  DEPENDS_ON mesh_from_inline_triangle mesh_from_inline_triangle_no_adapt
             mesh_from_inline_triangle_polygon
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving mesh_from_inline_triangle subdirectory")
# ------------------------------------------------------------------------------
