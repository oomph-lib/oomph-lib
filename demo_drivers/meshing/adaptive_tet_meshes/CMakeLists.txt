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
message(VERBOSE "Entered adaptive_tet_meshes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(adaptive_tet_meshes C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME geom_obj_with_boundary
  SOURCES geom_obj_with_boundary.cc
  LIBRARIES oomph::meshes oomph::generic)

if(OOMPH_HAS_CGAL)
  oomph_add_executable(
    NAME planar_facet_bounded_mesh_from_gmsh
    SOURCES planar_facet_bounded_tetmesh.cc tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson oomph::meshes oomph::generic)

  oomph_add_executable(
    NAME planar_facet_bounded_mesh_from_tetgen
    SOURCES planar_facet_bounded_tetmesh.cc tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson oomph::meshes oomph::generic
    CXX_DEFINITIONS DO_TETGEN)

  oomph_add_executable(
    NAME curved_facet_bounded_mesh_from_gmsh
    SOURCES curved_facet_bounded_tetmesh.cc tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson oomph::meshes oomph::generic)

  oomph_add_executable(
    NAME curved_facet_bounded_mesh_from_tetgen
    SOURCES curved_facet_bounded_tetmesh.cc tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson oomph::meshes oomph::generic
    CXX_DEFINITIONS DO_TETGEN)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS geom_obj_with_boundary)

if(OOMPH_HAS_CGAL)
  list(
    APPEND TEST_TARGET_DEPS planar_facet_bounded_mesh_from_gmsh
    planar_facet_bounded_mesh_from_tetgen curved_facet_bounded_mesh_from_gmsh
    curved_facet_bounded_mesh_from_tetgen)
endif()

oomph_add_test(
  TEST_NAME meshing.adaptive_tet_meshes
  DEPENDS_ON "${TEST_TARGET_DEPS}"
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_GMSH_COMMAND}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving adaptive_tet_meshes subdirectory")
# ------------------------------------------------------------------------------
