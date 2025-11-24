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
message(VERBOSE "Entered multi_physics subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(multi_physics C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    boussinesq_convection
    3d_boussinesq_convection
    double_diffusive_convection
    marangoni_convection
    insoluble_surfactant
    rayleigh_instability_surfactant
    soluble_surfactant
    two_layer_soluble_surfactant
    thermo
    thermal_fibre
    axisym_heat_sphere
    b_convection_sphere
    spherical_shell_convection)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving multi_physics subdirectory")
# ------------------------------------------------------------------------------
