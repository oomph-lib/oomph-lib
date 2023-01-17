# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_from_tetgen subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mesh_from_tetgen C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME mesh_from_tetgen_navier_stokes
  SOURCES mesh_from_tetgen_navier_stokes.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_tetgen_poisson
  SOURCES mesh_from_tetgen_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_tetgen
  TARGET_DEPENDENCIES mesh_from_tetgen_navier_stokes mesh_from_tetgen_poisson
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES cube_hole.1.ele cube_hole.1.face cube_hole.1.node
                 cube_hole.poly
  LABELS meshing mesh_from_tetgen)

message(VERBOSE "Leaving mesh_from_tetgen subdirectory")
# ------------------------------------------------------------------------------
