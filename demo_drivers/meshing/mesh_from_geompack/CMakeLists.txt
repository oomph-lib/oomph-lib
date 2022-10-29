# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_from_geompack subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(mesh_from_geompack C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME convert_geom_file
  SOURCES convert_geom_file.cc
  CXX_STANDARD 14
  SILENCE_NO_LIBS_SUPPLIED_WARNING)

oomph_add_executable(
  NAME mesh_from_geompack_poisson
  SOURCES mesh_from_geompack_poisson.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_geompack
  TARGET_DEPENDENCIES convert_geom_file mesh_from_geompack_poisson
  EXTRA_REQUIRES
    box_hole.cs2
    box_hole.m2
    box_hole.mh2
    box_hole.rg2
    quart_circ.cs2
    quart_circ.m2
    quart_circ.mh2
    quart_circ.rg2
    rect.cs2
    rect.m2
    rect.mh2
    rect.rg2
  LABELS meshing mesh_from_geompack)

message(VERBOSE "Leaving mesh_from_geompack subdirectory")
# ------------------------------------------------------------------------------
