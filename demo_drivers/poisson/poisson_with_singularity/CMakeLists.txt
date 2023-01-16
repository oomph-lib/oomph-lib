# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered poisson_with_singularity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(poisson_with_singularity C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson
  SOURCES two_d_poisson.cc poisson_elements_with_singularity.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson_without_singularity_subtracted_off
  SOURCES two_d_poisson.cc poisson_elements_with_singularity.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic
  CXX_DEFINITIONS DONT_USE_SINGULARITY)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.poisson_with_singularity
  TARGET_DEPENDENCIES two_d_poisson
                      two_d_poisson_without_singularity_subtracted_off
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS poisson poisson_with_singularity)

message(VERBOSE "Leaving poisson_with_singularity subdirectory")
# ------------------------------------------------------------------------------
