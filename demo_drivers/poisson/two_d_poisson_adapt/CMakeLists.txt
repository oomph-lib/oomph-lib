# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_poisson_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_poisson_adapt C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_adapt
  SOURCES two_d_poisson_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.two_d_poisson_adapt
  DEPENDS_ON two_d_poisson_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS poisson two_d_poisson_adapt adapt)

message(VERBOSE "Leaving two_d_poisson_adapt subdirectory")
# ------------------------------------------------------------------------------
