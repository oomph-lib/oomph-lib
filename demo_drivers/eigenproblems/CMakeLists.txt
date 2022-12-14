# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigenproblems subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eigenproblems C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS)

if(HAVE_ARPACK_SOURCES)
  list(APPEND SUBDIRS harmonic orr_sommerfeld)
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving eigenproblems subdirectory")
# ------------------------------------------------------------------------------
