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
message(VERBOSE "Entered meshing subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(meshing C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    mesh_from_inline_triangle
    mesh_from_inline_triangle_internal_boundaries
    mesh_from_triangle
    mesh_from_geompack
    mesh_from_tetgen
    mesh_from_xfig_triangle
    mesh_from_vmtk
    quad_from_triangle_mesh
    quadratic_vmtk_tetgen
    adaptive_tet_meshes)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving meshing subdirectory")
# ------------------------------------------------------------------------------
