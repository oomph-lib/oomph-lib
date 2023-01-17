# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered FAQ subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
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
  TARGET_DEPENDENCIES broken_this_demo customising_output
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS FAQ
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving FAQ subdirectory")
# ------------------------------------------------------------------------------
