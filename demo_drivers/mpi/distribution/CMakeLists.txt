# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered distribution subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(distribution C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
