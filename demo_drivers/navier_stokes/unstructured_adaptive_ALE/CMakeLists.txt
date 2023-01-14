# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_ALE subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_adaptive_ALE C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_adaptive_ALE
  SOURCES unstructured_adaptive_ALE.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::rigid_body oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_adaptive_ALE
  TARGET_DEPENDENCIES unstructured_adaptive_ALE
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS navier_stokes unstructured_adaptive_ALE adapt)

message(VERBOSE "Leaving unstructured_adaptive_ALE subdirectory")
# ------------------------------------------------------------------------------
