# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mesh_from_vmtk subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(mesh_from_vmtk C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
  SOURCES create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc
  LIBRARIES oomph::navier_stokes)

# Symlink built mesh to scripts directory
set(NAME create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh)
set(LINK ../../../../scripts/${TARGET})
add_custom_target(
  install_mesh
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${TARGET} ${LINK}
  DEPENDS ${TARGET})

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME meshing.mesh_from_vmtk
  TARGET_DEPENDENCIES create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
                      create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
  EXTRA_REQUIRES iliac.xda replies.txt
  LABELS meshing mesh_from_vmtk)

message(VERBOSE "Leaving mesh_from_vmtk subdirectory")
# ------------------------------------------------------------------------------
