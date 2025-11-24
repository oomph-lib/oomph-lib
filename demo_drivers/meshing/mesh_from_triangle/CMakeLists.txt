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
message(VERBOSE "Entered mesh_from_triangle subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(mesh_from_triangle C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME mesh_from_triangle_navier_stokes
  SOURCES mesh_from_triangle_navier_stokes.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_triangle_navier_stokes_adapt
  SOURCES mesh_from_triangle_navier_stokes.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPT)

oomph_add_executable(
  NAME mesh_from_triangle_poisson
  SOURCES mesh_from_triangle_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_triangle_poisson_adapt
  SOURCES mesh_from_triangle_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPT)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_triangle
  DEPENDS_ON mesh_from_triangle_navier_stokes
             mesh_from_triangle_navier_stokes_adapt mesh_from_triangle_poisson
             mesh_from_triangle_poisson_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             box_hole.1.edge
             box_hole.1.ele
             box_hole.1.node
             box_hole.1.poly
             box_hole.poly
             box_hole_adapt.1.ele
             box_hole_adapt.1.node
             box_hole_adapt.1.poly
             box_hole_adapt.poly
             flow_past_box.1.edge
             flow_past_box.1.ele
             flow_past_box.1.node
             flow_past_box.1.poly
             flow_past_box.poly
             flow_past_box_adapt.1.ele
             flow_past_box_adapt.1.node
             flow_past_box_adapt.1.poly
             flow_past_box_adapt.poly)

message(VERBOSE "Leaving mesh_from_triangle subdirectory")
# ------------------------------------------------------------------------------
