# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered oscillating_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(oscillating_shell C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME oscillating_shell
  SOURCES oscillating_shell.cc
  LIBRARIES oomph::shell oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.oscillating_shell
  DEPENDS_ON oscillating_shell
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving oscillating_shell subdirectory")
# ------------------------------------------------------------------------------
