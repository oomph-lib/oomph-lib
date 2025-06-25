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
message(VERBOSE "Entered mesh_from_geompack subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(mesh_from_geompack C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME convert_geom_file
  SOURCES convert_geom_file.cc
  SILENCE_NO_LIBS_SUPPLIED_WARNING)

oomph_add_executable(
  NAME mesh_from_geompack_poisson
  SOURCES mesh_from_geompack_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_geompack
  DEPENDS_ON convert_geom_file mesh_from_geompack_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
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
             rect.rg2)

message(VERBOSE "Leaving mesh_from_geompack subdirectory")
# ------------------------------------------------------------------------------
