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
message(VERBOSE "Entered spherical_shell_convection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(spherical_shell_convection C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME spherical_shell_convection
  SOURCES spherical_shell_convection.cc spherical_buoyant_navier_stokes.h
  LIBRARIES oomph::spherical_advection_diffusion oomph::spherical_navier_stokes
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.spherical_shell_convection
  DEPENDS_ON spherical_shell_convection
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving spherical_shell_convection subdirectory")
# ------------------------------------------------------------------------------
