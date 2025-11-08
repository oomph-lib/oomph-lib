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
message(VERBOSE "Entered double_diffusive_convection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(double_diffusive_convection C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME dd_convection
  SOURCES dd_convection.cc db_navier_st_elements.h
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME multimesh_dd_convection
  SOURCES multimesh_dd_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME multimesh_ref_dd_convection
  SOURCES multimesh_ref_dd_convection.cc db_nst_external_elements.h
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.double_diffusive_convection
  DEPENDS_ON dd_convection multimesh_dd_convection multimesh_ref_dd_convection
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving double_diffusive_convection subdirectory")
# ------------------------------------------------------------------------------
