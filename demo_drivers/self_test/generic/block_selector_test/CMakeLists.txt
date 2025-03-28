# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered block_selector_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(block_selector_test C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME block_selector_test
  SOURCES block_selector_test.cc
  LIBRARIES oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.generic.block_selector_test
  DEPENDS_ON block_selector_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
             compare_block_selector_nonnull_replacement_block_pt.sh)

message(VERBOSE "Leaving block_selector_test subdirectory")
# ------------------------------------------------------------------------------
