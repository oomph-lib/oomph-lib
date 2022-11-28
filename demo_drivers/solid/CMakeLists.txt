# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered solid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(solid C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
