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
message(VERBOSE "Entered axisym_navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(axisym_navier_stokes C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    axi_static_cap
    two_fluid_spherical_cap
    rotating_ends
    spin_up
    torus
    single_layer_free_surface_axisym
    two_layer_interface_axisym
    rayleigh_instability
    fibre
    counter_rotating_disks)

if(NOT OOMPH_SUPPRESSED_TRIANGLE_LIB)
  list(APPEND SUBDIRS unstructured_torus)
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
