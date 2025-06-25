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
message(VERBOSE "Entered locate_zeta subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(locate_zeta C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME locate_zeta_tester
  SOURCES locate_zeta_tester.cc locate_zeta_tester.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME locate_zeta_tester_3d
  SOURCES locate_zeta_tester_3d.cc locate_zeta_tester.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME locate_zeta_tester_triangle
  SOURCES locate_zeta_tester_triangle.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME locate_zeta_tester_tetgen
  SOURCES locate_zeta_tester_tetgen.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.locate_zeta
  DEPENDS_ON locate_zeta_tester locate_zeta_tester_3d
             locate_zeta_tester_triangle locate_zeta_tester_tetgen
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata cube_hole.1.node cube_hole.1.ele
             cube_hole.1.face)

message(VERBOSE "Leaving locate_zeta subdirectory")
# ------------------------------------------------------------------------------
