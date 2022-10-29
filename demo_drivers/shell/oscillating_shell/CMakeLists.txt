# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oscillating_shell subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(oscillating_shell C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME oscillating_shell
  SOURCES oscillating_shell.cc
  LIBRARIES oomph::shell
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.oscillating_shell
  TARGET_DEPENDENCIES oscillating_shell
  EXTRA_REQUIRES
  LABELS shell oscillating_shell)

message(VERBOSE "Leaving oscillating_shell subdirectory")
# ------------------------------------------------------------------------------
