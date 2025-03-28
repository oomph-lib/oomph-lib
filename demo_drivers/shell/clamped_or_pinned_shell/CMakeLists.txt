# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered clamped_or_pinned_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(clamped_or_pinned_shell C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME clamped_or_pinned_shell
  SOURCES clamped_or_pinned_shell.cc
  LIBRARIES oomph::shell oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.clamped_or_pinned_shell
  DEPENDS_ON clamped_or_pinned_shell
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving clamped_or_pinned_shell subdirectory")
# ------------------------------------------------------------------------------
