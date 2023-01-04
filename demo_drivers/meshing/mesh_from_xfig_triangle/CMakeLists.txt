# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_from_xfig_triangle subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mesh_from_xfig_triangle C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fig2poly
  SOURCES fig2poly.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_xfig_poisson
  SOURCES mesh_from_xfig_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME mesh_from_xfig_navier_stokes
  SOURCES mesh_from_xfig_navier_stokes.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# Create a symlink to the built mesh from the scripts directory
#
# TODO: Ask whether it's better to make this a copy rather a symlink since the
# build/ directory is likely to be deleted
add_custom_target(
  install_fig2poly COMMAND ln -sf "${CMAKE_BINARY_DIR}/fig2poly"
                           "${OOMPH_ROOT_DIR}/scripts/fig2poly")
oomph_get_hashed_target_name(fig2poly FIG2POLY_HASHED_TARGET_NAME)
add_dependencies(install_fig2poly ${FIG2POLY_HASHED_TARGET_NAME})

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_xfig_triangle
  TARGET_DEPENDENCIES fig2poly mesh_from_xfig_poisson
                      mesh_from_xfig_navier_stokes
  EXTRA_REQUIRES
    hole.fig
    hole.fig.1.ele
    hole.fig.1.node
    hole.fig.1.poly
    hole.fig.poly
    oomph_mesh.fig
    oomph_mesh.fig.1.ele
    oomph_mesh.fig.1.node
    oomph_mesh.fig.1.poly
    oomph_mesh.fig.poly
  LABELS meshing mesh_from_xfig_triangle)

message(VERBOSE "Leaving mesh_from_xfig_triangle subdirectory")
# ------------------------------------------------------------------------------
