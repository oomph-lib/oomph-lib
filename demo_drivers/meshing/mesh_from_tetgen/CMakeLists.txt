# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_from_tetgen subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(mesh_from_tetgen C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME mesh_from_tetgen_navier_stokes
  SOURCES mesh_from_tetgen_navier_stokes.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME mesh_from_tetgen_poisson
  SOURCES mesh_from_tetgen_poisson.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_tetgen
  TARGET_DEPENDENCIES mesh_from_tetgen_navier_stokes mesh_from_tetgen_poisson
  EXTRA_REQUIRES cube_hole.1.ele cube_hole.1.face cube_hole.1.node
                 cube_hole.poly
  LABELS "meshing;mesh_from_tetgen")

message(VERBOSE "Leaving mesh_from_tetgen subdirectory")
# ------------------------------------------------------------------------------
