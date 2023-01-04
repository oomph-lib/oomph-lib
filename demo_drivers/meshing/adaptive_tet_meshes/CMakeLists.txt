# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_tet_meshes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(adaptive_tet_meshes C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
  TARGET_DEPENDENCIES "${TEST_TARGET_DEPS}"
  EXTRA_REQUIRES
  LABELS meshing adaptive_tet_meshes adapt cgal)

message(VERBOSE "Leaving adaptive_tet_meshes subdirectory")
# ------------------------------------------------------------------------------
