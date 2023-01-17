# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multi_domain subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(multi_domain C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

set(SUBDIRS
    boussinesq_convection
    fsi_collapsible_channel
    fsi_osc_ring
    fsi_channel_with_leaflet
    turek_flag
    pseudo_solid_collapsible_tube)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving multi_domain subdirectory")
# ------------------------------------------------------------------------------
