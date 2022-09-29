# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_tet_meshes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(adaptive_tet_meshes C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME geom_obj_with_boundary
  SOURCES geom_obj_with_boundary.cc
  LIBRARIES oomph::generic
  CXX_STANDARD 14)

if(OOMPH_HAS_CGAL)
  oomph_add_executable(
    NAME planar_facet_bounded_mesh_from_gmsh
    SOURCES simple_shear.cc tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson
    CXX_STANDARD 14)

  oomph_add_executable(
    NAME planar_facet_bounded_mesh_from_tetgen
    SOURCES refineable_simple_shear.cc planar_facet_bounded_tetmesh.cc
            tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson
    CXX_DEFINITIONS DO_TETGEN
    CXX_STANDARD 14)

  oomph_add_executable(
    NAME curved_facet_bounded_mesh_from_gmsh
    SOURCES refineable_simple_shear.cc tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson
    CXX_STANDARD 14)

  oomph_add_executable(
    NAME curved_facet_bounded_mesh_from_tetgen
    SOURCES refineable_simple_shear.cc curved_facet_bounded_tetmesh.cc
            tetmesh_faceted_surfaces.h
    LIBRARIES oomph::poisson
    CXX_DEFINITIONS DO_TETGEN
    CXX_STANDARD 14)
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
  TARGET_DEPENDENCIES ${TEST_TARGET_DEPS}
  EXTRA_REQUIRES
  LABELS "meshing;adaptive_tet_meshes;adapt")

message(VERBOSE "Leaving adaptive_tet_meshes subdirectory")
# ------------------------------------------------------------------------------
