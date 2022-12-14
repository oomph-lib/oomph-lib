# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered locate_zeta subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(locate_zeta C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
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
  TEST_NAME locate_zeta
  TARGET_DEPENDENCIES locate_zeta_tester locate_zeta_tester_3d
                      locate_zeta_tester_triangle locate_zeta_tester_tetgen
  EXTRA_REQUIRES cube_hole.1.node cube_hole.1.ele cube_hole.1.face
  LABELS locate_zeta)

message(VERBOSE "Leaving locate_zeta subdirectory")
# ------------------------------------------------------------------------------
