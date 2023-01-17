# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered stored_shape_fcts subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(stored_shape_fcts C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_stored_shape_fcts
  SOURCES two_d_poisson_stored_shape_fcts.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.stored_shape_fcts
  TARGET_DEPENDENCIES two_d_poisson_stored_shape_fcts
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS optimisation stored_shape_fcts)

message(VERBOSE "Leaving stored_shape_fcts subdirectory")
# ------------------------------------------------------------------------------
