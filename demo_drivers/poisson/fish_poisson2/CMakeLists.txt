# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson2 subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fish_poisson2 C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fish_domain
  SOURCES fish_domain.cc
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fish_poisson_adapt
  SOURCES fish_poisson_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fish_poisson_no_adapt
  SOURCES fish_poisson_no_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fish_poisson_node_update
  SOURCES fish_poisson_node_update.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fish_poisson_simple_adapt
  SOURCES fish_poisson_simple_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.fish_poisson2
  DEPENDS_ON fish_domain fish_poisson_adapt fish_poisson_no_adapt
             fish_poisson_node_update fish_poisson_simple_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving fish_poisson2 subdirectory")
# ------------------------------------------------------------------------------
