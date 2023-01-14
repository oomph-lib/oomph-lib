# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_poisson_flux_bc_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_poisson_flux_bc_adapt C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_flux_bc_adapt
  SOURCES two_d_poisson_flux_bc_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.two_d_poisson_flux_bc_adapt
  TARGET_DEPENDENCIES two_d_poisson_flux_bc_adapt
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS poisson two_d_poisson_flux_bc_adapt adapt)

message(VERBOSE "Leaving two_d_poisson_flux_bc_adapt subdirectory")
# ------------------------------------------------------------------------------
