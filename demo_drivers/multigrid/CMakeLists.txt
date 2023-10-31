# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multigrid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(multigrid C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS scalar_multigrid helmholtz_multigrid)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving multigrid subdirectory")
# ------------------------------------------------------------------------------
