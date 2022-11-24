# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered clamped_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(clamped_shell C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME clamped_shell
  SOURCES clamped_shell.cc
  LIBRARIES oomph::shell)

oomph_add_executable(
  NAME clamped_shell_with_arclength_cont
  SOURCES clamped_shell_with_arclength_cont.cc
  LIBRARIES oomph::shell)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.clamped_shell
  TARGET_DEPENDENCIES clamped_shell clamped_shell_with_arclength_cont
  EXTRA_REQUIRES
  LABELS shell clamped_shell)

message(VERBOSE "Leaving clamped_shell subdirectory")
# ------------------------------------------------------------------------------
