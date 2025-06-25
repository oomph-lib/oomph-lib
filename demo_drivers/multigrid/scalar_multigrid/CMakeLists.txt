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
message(VERBOSE "Entered scalar_multigrid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(scalar_multigrid C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME two_d_poisson_tanh_flux_bc
  SOURCES two_d_poisson_tanh_flux_bc.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson_tanh_flux_bc_validate
  SOURCES two_d_poisson_tanh_flux_bc_validate.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unit_cube_poisson
  SOURCES unit_cube_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unit_cube_poisson_validate
  SOURCES unit_cube_poisson_validate.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multigrid.scalar_multigrid
  DEPENDS_ON two_d_poisson_tanh_flux_bc two_d_poisson_tanh_flux_bc_validate
             unit_cube_poisson unit_cube_poisson_validate
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving scalar_multigrid subdirectory")
# ------------------------------------------------------------------------------
