# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered free_surface_rotation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(free_surface_rotation C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME free_surface_rotation
  SOURCES free_surface_rotation.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.free_surface_rotation
  DEPENDS_ON free_surface_rotation
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving free_surface_rotation subdirectory")
# ------------------------------------------------------------------------------
