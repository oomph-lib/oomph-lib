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
message(VERBOSE "Entered quadratic_vmtk_tetgen subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(quadratic_vmtk_tetgen C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME snap_mesh
  SOURCES snap_mesh.cc
  LIBRARIES oomph::mesh_smoothing
            oomph::solid
            oomph::constitutive
            oomph::navier_stokes
            oomph::poisson
            oomph::linear_elasticity
            oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME xda_to_poly_fsi
  SOURCES xda_to_poly_fsi.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.quadratic_vmtk_tetgen
  DEPENDS_ON snap_mesh xda_to_poly_fsi
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata short_coarse_mesh_files)

# This test needs more than the default allotted time to complete
set_tests_properties(meshing.quadratic_vmtk_tetgen PROPERTIES TIMEOUT 3600)

message(VERBOSE "Leaving quadratic_vmtk_tetgen subdirectory")
# ------------------------------------------------------------------------------
