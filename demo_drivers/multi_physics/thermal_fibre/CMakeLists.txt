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
message(VERBOSE "Entered thermal_fibre subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(thermal_fibre C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME thermal_fibre
  SOURCES non_isothermal_axisym_nst_elements.h thermal_fibre.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface
            oomph::steady_axisym_advection_diffusion oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.thermal_fibre
  DEPENDS_ON thermal_fibre
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving thermal_fibre subdirectory")
# ------------------------------------------------------------------------------
