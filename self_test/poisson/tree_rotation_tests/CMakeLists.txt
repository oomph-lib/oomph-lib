# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered tree_rotation_tests subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(tree_rotation_tests C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME tree_2d
  SOURCES tree_2d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME tree_3d
  SOURCES tree_3d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME tree_rotation_tests
  TARGET_DEPENDENCIES tree_2d tree_3d
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS poisson tree_rotation
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving tree_rotation_tests subdirectory")
# ------------------------------------------------------------------------------
