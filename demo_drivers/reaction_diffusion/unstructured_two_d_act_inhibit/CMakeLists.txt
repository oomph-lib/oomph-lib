# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_two_d_act_inhibit subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_two_d_act_inhibit C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_act_inhibit
  SOURCES unstructured_two_d_act_inhibit.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME reaction_diffusion.unstructured_two_d_act_inhibit
  DEPENDS_ON unstructured_two_d_act_inhibit
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unstructured_two_d_act_inhibit subdirectory")
# ------------------------------------------------------------------------------
