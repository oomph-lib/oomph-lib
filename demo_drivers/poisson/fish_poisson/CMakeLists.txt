# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fish_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fish_poisson
  SOURCES fish_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fish_poisson_no_adapt
  SOURCES fish_poisson_no_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.fish_poisson
  DEPENDS_ON fish_poisson fish_poisson_no_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving fish_poisson subdirectory")
# ------------------------------------------------------------------------------
