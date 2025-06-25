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
message(VERBOSE "Entered mesh_from_vmtk subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(mesh_from_vmtk C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
  SOURCES create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

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
  DEPENDS_ON create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
             create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata iliac.xda replies.txt)

message(VERBOSE "Leaving mesh_from_vmtk subdirectory")
# ------------------------------------------------------------------------------
