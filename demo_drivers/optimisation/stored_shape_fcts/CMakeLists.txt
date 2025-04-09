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
  DEPENDS_ON two_d_poisson_stored_shape_fcts
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving stored_shape_fcts subdirectory")
# ------------------------------------------------------------------------------
