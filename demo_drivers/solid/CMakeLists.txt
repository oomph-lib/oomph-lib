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
message(VERBOSE "Entered solid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(solid C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    disk_compression
    disk_oscillation
    shock_disk
    static_fish
    airy_cantilever
    three_d_cantilever
    simple_shear
    compressed_square
    prescribed_displ_lagr_mult
    three_d_prescribed_displ_lagr_mult
    three_d_traction
    pseudo_solid_collapsible_channel
    vmtk_solid
    unstructured_solid
    unstructured_three_d_solid
    unstructured_adaptive_solid)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving solid subdirectory")
# ------------------------------------------------------------------------------
