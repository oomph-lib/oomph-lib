# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_fourier_decomposed_helmholtz subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(pml_fourier_decomposed_helmholtz C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS oscillating_sphere)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving pml_fourier_decomposed_helmholtz subdirectory")
# ------------------------------------------------------------------------------
