# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flux_control subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(flux_control C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME flux_control
  SOURCES flux_control.cc
  LIBRARIES oomph::womersley oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.flux_control
  DEPENDS_ON flux_control
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving flux_control subdirectory")
# ------------------------------------------------------------------------------
