# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(poisson C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS convergence_tests tree_rotation_tests octree_test face_tests)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving poisson subdirectory")
# ------------------------------------------------------------------------------
