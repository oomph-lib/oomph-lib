# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered heat_transfer_and_melting subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(heat_transfer_and_melting C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS two_d_unsteady_heat_melt)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving heat_transfer_and_melting subdirectory")
# ------------------------------------------------------------------------------
