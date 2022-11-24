# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_bretherton subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_bretherton C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_breth
  SOURCES three_d_breth.cc
          canyon_spine_mesh.h
          comb_tip_spine_mesh.h
          merge_meshes.h
          tip_spine_mesh.h
          comb_can_spine_mesh.h
          extra_elements.h
          st_mesh.h
  LIBRARIES oomph::navier_stokes oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.three_d_bretherton
  TARGET_DEPENDENCIES three_d_breth
  EXTRA_REQUIRES
  LABELS navier_stokes three_d_bretherton)

message(VERBOSE "Leaving three_d_bretherton subdirectory")
# ------------------------------------------------------------------------------
