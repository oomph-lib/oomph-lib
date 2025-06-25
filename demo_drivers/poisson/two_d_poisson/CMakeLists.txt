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
message(VERBOSE "Entered two_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_poisson C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME two_d_poisson
  SOURCES two_d_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson2
  SOURCES two_d_poisson2.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# Do we have the hsl sources?
set(CXX_DEFINITION)
if(HAVE_HSL_SOURCES)
  set(CXX_DEFINITION HAVE_HSL_SOURCES)
endif()

oomph_add_executable(
  NAME two_d_poisson_compare_solvers
  SOURCES two_d_poisson_compare_solvers.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic
  CXX_DEFINITIONS ${CXX_DEFINITION})

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.two_d_poisson
  DEPENDS_ON two_d_poisson two_d_poisson2 two_d_poisson_compare_solvers
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving two_d_poisson subdirectory")
# ------------------------------------------------------------------------------
