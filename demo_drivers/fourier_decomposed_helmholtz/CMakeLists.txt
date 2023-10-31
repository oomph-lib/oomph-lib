# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fourier_decomposed_helmholtz subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fourier_decomposed_helmholtz C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS sphere_scattering)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving fourier_decomposed_helmholtz subdirectory")
# ------------------------------------------------------------------------------
