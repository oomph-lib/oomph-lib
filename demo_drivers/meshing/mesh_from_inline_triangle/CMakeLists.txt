# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_from_inline_triangle subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mesh_from_inline_triangle C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
  TARGET_DEPENDENCIES
    mesh_from_inline_triangle mesh_from_inline_triangle_no_adapt
    mesh_from_inline_triangle_polygon
  EXTRA_REQUIRES
  LABELS meshing mesh_from_inline_triangle)

message(VERBOSE "Leaving mesh_from_inline_triangle subdirectory")
# ------------------------------------------------------------------------------
