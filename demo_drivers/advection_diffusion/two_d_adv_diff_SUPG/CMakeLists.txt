# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_adv_diff_SUPG subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_adv_diff_SUPG C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_adv_diff_SUPG
  SOURCES two_d_adv_diff_SUPG.cc
  LIBRARIES oomph::advection_diffusion oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME advection_diffusion.two_d_adv_diff_SUPG
  TARGET_DEPENDENCIES two_d_adv_diff_SUPG
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS advection_diffusion two_d_adv_diff_SUPG)

message(VERBOSE "Leaving two_d_adv_diff_SUPG subdirectory")
# ------------------------------------------------------------------------------
