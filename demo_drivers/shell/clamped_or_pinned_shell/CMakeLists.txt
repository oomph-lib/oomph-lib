# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered clamped_or_pinned_shell subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(clamped_or_pinned_shell C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME clamped_or_pinned_shell
  SOURCES clamped_or_pinned_shell.cc
  LIBRARIES oomph::shell
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.clamped_or_pinned_shell
  TARGET_DEPENDENCIES clamped_or_pinned_shell
  EXTRA_REQUIRES
  LABELS shell clamped_or_pinned_shell)

message(VERBOSE "Leaving clamped_or_pinned_shell subdirectory")
# ------------------------------------------------------------------------------
