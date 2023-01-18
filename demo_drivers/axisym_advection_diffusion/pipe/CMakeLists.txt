# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pipe subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(pipe C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME pipe
  SOURCES pipe.cc
  LIBRARIES oomph::axisym_advection_diffusion oomph::meshes oomph::generic)

oomph_add_executable(
  NAME pipe_variable_diff
  SOURCES pipe_variable_diff.cc
  LIBRARIES oomph::axisym_advection_diffusion oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_advection_diffusion.pipe
  DEPENDS_ON pipe pipe_variable_diff
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS axisym_advection_diffusion pipe)

message(VERBOSE "Leaving pipe subdirectory")
# ------------------------------------------------------------------------------
