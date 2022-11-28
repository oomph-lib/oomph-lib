# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axisym_navier_stokes C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
