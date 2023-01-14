# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_poisson C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson
  SOURCES two_d_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson2
  SOURCES two_d_poisson2.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# Do we have the hsl sources?
set(CXX_DEFINITION)
if(HAVE_HSL_SOURCES)
  set(CXX_DEFINITION HAVE_HSL_SOURCES)
endif()

oomph_add_executable(
  NAME two_d_poisson_compare_solvers
  SOURCES two_d_poisson_compare_solvers.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic
  CXX_DEFINITIONS ${CXX_DEFINITION})

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.two_d_poisson
  TARGET_DEPENDENCIES two_d_poisson two_d_poisson2 two_d_poisson_compare_solvers
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS poisson two_d_poisson)

message(VERBOSE "Leaving two_d_poisson subdirectory")
# ------------------------------------------------------------------------------
