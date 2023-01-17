# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(shell C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS clamped_shell clamped_or_pinned_shell oscillating_shell plate)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving shell subdirectory")
# ------------------------------------------------------------------------------
