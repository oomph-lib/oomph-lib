# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered advection_diffusion subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(advection_diffusion C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS two_d_adv_diff_flux_bc two_d_adv_diff_adapt two_d_adv_diff_SUPG
    two_d_variable_diff)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
