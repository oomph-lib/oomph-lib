# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered faceted_surface_vtu_output subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(faceted_surface_vtu_output C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME faceted_surface_vtu_output_tester
  SOURCES faceted_surface_vtu_output_tester.cc
  LIBRARIES oomph::poisson oomph::generic)

message(VERBOSE "Leaving faceted_surface_vtu_output subdirectory")
# ------------------------------------------------------------------------------
