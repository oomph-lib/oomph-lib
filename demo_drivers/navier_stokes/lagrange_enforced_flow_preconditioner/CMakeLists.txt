# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered lagrange_enforced_flow_preconditioner subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(lagrange_enforced_flow_preconditioner C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_tilted_square
  SOURCES two_d_tilted_square.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.lagrange_enforced_flow_preconditioner
  DEPENDS_ON two_d_tilted_square
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes lagrange_enforced_flow_preconditioner)

message(VERBOSE "Leaving lagrange_enforced_flow_preconditioner subdirectory")
# ------------------------------------------------------------------------------
