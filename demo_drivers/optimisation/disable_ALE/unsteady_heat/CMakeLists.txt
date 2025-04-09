# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_heat subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unsteady_heat C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_unsteady_heat
  SOURCES two_d_unsteady_heat.cc
  LIBRARIES oomph::unsteady_heat oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_unsteady_heat_adapt
  SOURCES two_d_unsteady_heat_adapt.cc
  LIBRARIES oomph::unsteady_heat oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.unsteady_heat
  DEPENDS_ON two_d_unsteady_heat two_d_unsteady_heat_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
