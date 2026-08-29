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
message(VERBOSE "Entered shell subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(shell C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME shell
  SOURCES shell.cc benchmark_elements.cc
  LIBRARIES oomph::axisym_solid oomph::constitutive oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME shell_adaptive
  SOURCES shell_adaptive.cc benchmark_elements.cc
  LIBRARIES oomph::axisym_solid oomph::constitutive oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_solid.shell
  DEPENDS_ON shell shell_adaptive
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving shell subdirectory")
# ------------------------------------------------------------------------------
