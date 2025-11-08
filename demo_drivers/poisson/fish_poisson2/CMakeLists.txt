# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson2 subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
