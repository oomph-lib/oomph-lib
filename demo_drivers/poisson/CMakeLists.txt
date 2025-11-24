# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(poisson C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    poisson_with_singularity
    one_d_poisson
    one_d_poisson_generic_only
    one_d_poisson_adapt
    one_d_poisson_hp_adapt
    two_d_poisson
    two_d_poisson_flux_bc
    two_d_poisson_flux_bc2
    two_d_poisson_flux_bc_adapt
    fish_poisson
    two_d_poisson_adapt
    two_d_poisson_hp_adapt
    fish_poisson2
    elastic_poisson
    eighth_sphere_poisson
    eighth_sphere_poisson_hp_adapt
    spectral)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving poisson subdirectory")
# ------------------------------------------------------------------------------
