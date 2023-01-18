# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered clamped_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(clamped_shell C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME clamped_shell
  SOURCES clamped_shell.cc
  LIBRARIES oomph::shell oomph::meshes oomph::generic)

oomph_add_executable(
  NAME clamped_shell_with_arclength_cont
  SOURCES clamped_shell_with_arclength_cont.cc
  LIBRARIES oomph::shell oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME shell.clamped_shell
  DEPENDS_ON clamped_shell clamped_shell_with_arclength_cont
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS shell clamped_shell)

message(VERBOSE "Leaving clamped_shell subdirectory")
# ------------------------------------------------------------------------------
