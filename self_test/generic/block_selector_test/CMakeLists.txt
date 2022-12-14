# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered block_selector_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(block_selector_test C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME block_selector_test
  SOURCES block_selector_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME generic.block_selector_test
  TARGET_DEPENDENCIES block_selector_test
  EXTRA_REQUIRES compare_block_selector_nonnull_replacement_block_pt.sh
  LABELS block_selector_test)

message(VERBOSE "Leaving block_selector_test subdirectory")
# ------------------------------------------------------------------------------
