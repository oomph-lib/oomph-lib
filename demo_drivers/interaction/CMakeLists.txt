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
message(VERBOSE "Entered interaction subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(interaction C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    free_boundary_poisson
    fsi_osc_ring
    fsi_collapsible_channel
    fsi_channel_with_leaflet
    turek_flag
    fsi_driven_cavity
    fsi_channel_seg_and_precond
    vmtk_fsi
    unsteady_vmtk_fsi
    unstructured_fsi
    unstructured_three_d_fsi
    unstructured_adaptive_fsi
    elastic_bretherton
    pseudo_solid_collapsible_tube
    linearised_fsi_pulsewave
    linearised_poroelastic_fsi_pulsewave
    pseudo_solid_fsi_channel_with_leaflet
    acoustic_fsi
    fourier_decomposed_acoustic_fsi)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving interaction subdirectory")
# ------------------------------------------------------------------------------
