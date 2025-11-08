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
message(VERBOSE "Entered FAQ subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(FAQ C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME broken_this_demo
  SOURCES broken_this_demo.cc
  SILENCE_NO_LIBS_SUPPLIED_WARNING)

oomph_add_executable(
  NAME customising_output
  SOURCES customising_output.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME FAQ.FAQ
  DEPENDS_ON broken_this_demo customising_output
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving FAQ subdirectory")
# ------------------------------------------------------------------------------
