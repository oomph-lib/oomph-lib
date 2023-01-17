# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered insoluble_surfactant subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(insoluble_surfactant C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME surfactant
  SOURCES surfactant.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
            oomph::multi_physics oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.insoluble_surfactant
  TARGET_DEPENDENCIES surfactant
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS multi_physics insoluble_surfactant)

message(VERBOSE "Leaving insoluble_surfactant subdirectory")
# ------------------------------------------------------------------------------
