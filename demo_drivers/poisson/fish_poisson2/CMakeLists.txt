# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson2 subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(fish_poisson2 C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME fish_domain
  SOURCES fish_domain.cc
  LIBRARIES oomph::generic
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fish_poisson_adapt
  SOURCES fish_poisson_adapt.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fish_poisson_no_adapt
  SOURCES fish_poisson_no_adapt.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fish_poisson_node_update
  SOURCES fish_poisson_node_update.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME fish_poisson_simple_adapt
  SOURCES fish_poisson_simple_adapt.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.fish_poisson2
  TARGET_DEPENDENCIES fish_domain fish_poisson_adapt fish_poisson_no_adapt
                      fish_poisson_node_update fish_poisson_simple_adapt
  EXTRA_REQUIRES
  LABELS "poisson;fish_poisson2")

message(VERBOSE "Leaving fish_poisson2 subdirectory")
# ------------------------------------------------------------------------------
