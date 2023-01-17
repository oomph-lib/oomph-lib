# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered reaction_diffusion subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(reaction_diffusion C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS one_d_act_inhibit two_d_act_inhibit two_d_act_inhibit_adapt
    unstructured_two_d_act_inhibit)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving reaction_diffusion subdirectory")
# ------------------------------------------------------------------------------
