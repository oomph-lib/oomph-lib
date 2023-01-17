# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eighth_sphere_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eighth_sphere_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME eighth_sphere_poisson
  SOURCES eighth_sphere_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.eighth_sphere_poisson
  TARGET_DEPENDENCIES eighth_sphere_poisson
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS poisson eighth_sphere_poisson)

message(VERBOSE "Leaving eighth_sphere_poisson subdirectory")
# ------------------------------------------------------------------------------
