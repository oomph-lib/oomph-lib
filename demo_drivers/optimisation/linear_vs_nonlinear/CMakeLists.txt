# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linear_vs_nonlinear subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(linear_vs_nonlinear C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson
  SOURCES two_d_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.linear_vs_nonlinear
  TARGET_DEPENDENCIES two_d_poisson
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS optimisation linear_vs_nonlinear)

message(VERBOSE "Leaving linear_vs_nonlinear subdirectory")
# ------------------------------------------------------------------------------
