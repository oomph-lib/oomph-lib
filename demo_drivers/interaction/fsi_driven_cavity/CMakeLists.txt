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
message(VERBOSE "Entered fsi_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_driven_cavity C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME fsi_driven_cavity_driver
  SOURCES fsi_driven_cavity_driver.cc fsi_driven_cavity_problem.h
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_driven_cavity_segregated_driver
  SOURCES fsi_driven_cavity_segregated_driver.cc fsi_driven_cavity_problem.h
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_driven_cavity
  DEPENDS_ON fsi_driven_cavity_driver fsi_driven_cavity_segregated_driver
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving fsi_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
