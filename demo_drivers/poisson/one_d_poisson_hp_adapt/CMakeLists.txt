# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_poisson_hp_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(one_d_poisson_hp_adapt C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_poisson_hp_adapt
  SOURCES one_d_poisson_hp_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.one_d_poisson_hp_adapt
  DEPENDS_ON one_d_poisson_hp_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS poisson one_d_poisson_hp_adapt adapt)

message(VERBOSE "Leaving one_d_poisson_hp_adapt subdirectory")
# ------------------------------------------------------------------------------
