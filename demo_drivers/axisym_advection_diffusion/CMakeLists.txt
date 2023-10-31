# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_advection_diffusion subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axisym_advection_diffusion C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS pipe)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
