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
message(VERBOSE "Entered mesh_from_tetgen subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
  DEPENDS_ON mesh_from_tetgen_navier_stokes mesh_from_tetgen_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata cube_hole.1.ele cube_hole.1.face
             cube_hole.1.node cube_hole.poly)

message(VERBOSE "Leaving mesh_from_tetgen subdirectory")
# ------------------------------------------------------------------------------
