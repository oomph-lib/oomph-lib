# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered turek_flag_non_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(turek_flag_non_fsi C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME turek_flag_non_fsi_alg
  SOURCES turek_flag_non_fsi.cc
  LIBRARIES oomph::navier_stokes)

oomph_add_executable(
  NAME turek_flag_non_fsi_macro
  SOURCES turek_flag_non_fsi.cc
  LIBRARIES oomph::navier_stokes
  CXX_DEFINITIONS USE_MACRO_ELEMENTS)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.turek_flag_non_fsi
  TARGET_DEPENDENCIES turek_flag_non_fsi_alg turek_flag_non_fsi_macro
  EXTRA_REQUIRES
  LABELS navier_stokes turek_flag_non_fsi)

message(VERBOSE "Leaving turek_flag_non_fsi subdirectory")
# ------------------------------------------------------------------------------
