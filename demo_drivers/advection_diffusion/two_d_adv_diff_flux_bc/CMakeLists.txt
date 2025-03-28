# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_adv_diff_flux_bc subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_adv_diff_flux_bc C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_adv_diff_flux_bc
  SOURCES two_d_adv_diff_flux_bc.cc
  LIBRARIES oomph::advection_diffusion oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME advection_diffusion.two_d_adv_diff_flux_bc
  DEPENDS_ON two_d_adv_diff_flux_bc
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving two_d_adv_diff_flux_bc subdirectory")
# ------------------------------------------------------------------------------
