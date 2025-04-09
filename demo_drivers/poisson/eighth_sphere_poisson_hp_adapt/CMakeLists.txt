# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eighth_sphere_poisson_hp_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eighth_sphere_poisson_hp_adapt C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME eighth_sphere_poisson_hp_adapt
  SOURCES eighth_sphere_poisson_hp_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.eighth_sphere_poisson_hp_adapt
  DEPENDS_ON eighth_sphere_poisson_hp_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving eighth_sphere_poisson_hp_adapt subdirectory")
# ------------------------------------------------------------------------------
