# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_advection_diffusion subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spherical_advection_diffusion C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS eluting_sphere)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving spherical_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
