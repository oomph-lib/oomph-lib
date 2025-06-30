# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered scripts subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(scripts C CXX Fortran)
if(PROJECT_IS_TOP_LEVEL)
  find_package(oomphlib CONFIG REQUIRED PATHS "../install")
endif()

# Build fig2poly from demo_drivers/meshing/mesh_from_xfig_triangle and get it to
# install the command
add_custom_target(
  fig2poly ALL
  COMMAND ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR} -B build
  COMMAND cd build
  COMMAND ${CMAKE_MAKE_PROGRAM} install_fig2poly
  WORKING_DIRECTORY
    "${CMAKE_CURRENT_LIST_DIR}/../demo_drivers/meshing/mesh_from_xfig_triangle")

message(VERBOSE "Leaving scripts subdirectory")
# ------------------------------------------------------------------------------
