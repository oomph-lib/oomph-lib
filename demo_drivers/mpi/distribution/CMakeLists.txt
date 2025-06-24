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
message(VERBOSE "Entered distribution subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(distribution C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

set(SUBDIRS
    hanging_node_reconciliation
    restart
    adaptive_driven_cavity
    hp_adaptive_driven_cavity
    circular_driven_cavity
    fish_poisson
    airy_cantilever
    two_d_poisson_flux_bc_adapt
    two_d_unstructured_adaptive_poisson
    hp_adaptive_poisson
    prescribed_displ_lagr_mult
    three_d_entry_flow
    three_d_cantilever
    three_d_prescribed_displ_lagr_mult
    clamped_shell
    bifurcation_tracking
    eigenproblem)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving distribution subdirectory")
# ------------------------------------------------------------------------------
