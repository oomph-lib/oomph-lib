# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered space_time_flow_past_oscillating_cylinder subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(space_time_oscillating_cylinder C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME space_time_oscillating_cylinder
  SOURCES space_time_oscillating_cylinder.cc
  LIBRARIES oomph::space_time_navier_stokes
            oomph::space_time_block_preconditioner oomph::navier_stokes
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.space_time_oscillating_cylinder
  DEPENDS_ON space_time_oscillating_cylinder
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE
        "Leaving space_time_flow_past_oscillating_cylinder subdirectory")
# ------------------------------------------------------------------------------
