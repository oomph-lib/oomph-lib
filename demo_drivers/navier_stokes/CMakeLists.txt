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
message(VERBOSE "Entered navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(navier_stokes C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    adaptive_driven_cavity
    adaptive_interface
    bretherton
    channel_with_leaflet
    circular_driven_cavity
    collapsible_channel
    curved_pipe
    driven_cavity
    falling_jet
    flow_past_cylinder
    flow_past_oscillating_cylinder
    flux_control
    free_surface_rotation
    hp_adaptive_driven_cavity
    inclined_plane
    lagrange_enforced_flow_preconditioner
    navier_stokes_with_singularity
    osc_ellipse
    osc_ring
    rayleigh_channel
    rayleigh_traction_channel
    schur_complement_preconditioner
    single_layer_free_surface
    space_time_flow_past_oscillating_cylinder
    spine_channel
    static_cap
    three_d_bretherton
    three_d_entry_flow
    three_d_free_surface
    three_d_static_cap
    turek_flag_non_fsi
    two_layer_interface
    unstructured_fluid
    unstructured_three_d_fluid
    vmtk_fluid
    vorticity_smoother)

if(NOT OOMPH_SUPPRESSED_TRIANGLE_LIB)
  list(APPEND SUBDIRS unstructured_adaptive_ALE unstructured_adaptive_fs
       jeffery_orbit)
endif()

if(NOT OOMPH_SUPPRESSED_TETGEN_LIB)
  list(APPEND SUBDIRS unstructured_adaptive_3d_ALE)
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving navier_stokes subdirectory")
# ------------------------------------------------------------------------------
