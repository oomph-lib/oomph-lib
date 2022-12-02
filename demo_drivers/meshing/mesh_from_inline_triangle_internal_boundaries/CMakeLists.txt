# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered mesh_from_inline_triangle_internal_boundaries subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mesh_from_inline_triangle_internal_boundaries C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME mesh_from_inline_triangle_internal_boundaries
  SOURCES mesh_from_inline_triangle_internal_boundaries.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_inline_triangle_internal_boundaries_extra
  SOURCES mesh_from_inline_triangle_internal_boundaries_extra.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_inline_triangle_internal_boundaries
  TARGET_DEPENDENCIES mesh_from_inline_triangle_internal_boundaries
                      mesh_from_inline_triangle_internal_boundaries_extra
  EXTRA_REQUIRES
  LABELS meshing mesh_from_inline_triangle_internal_boundaries)

message(VERBOSE
        "Leaving mesh_from_inline_triangle_internal_boundaries subdirectory")
# ------------------------------------------------------------------------------
