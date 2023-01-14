# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_3d_ALE subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_adaptive_3d_ALE C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

if(OOMPH_HAS_CGAL)
  oomph_add_executable(
    NAME uns_adapt_3d
    SOURCES uns_adapt_3d.cc
    LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
              oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_CGAL)
  oomph_add_test(
    TEST_NAME navier_stokes.unstructured_adaptive_3d_ALE
    TARGET_DEPENDENCIES uns_adapt_3d
    VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
    EXTRA_REQUIRES
    LABELS navier_stokes unstructured_adaptive_3d_ALE adapt)
endif()

message(VERBOSE "Leaving unstructured_adaptive_3d_ALE subdirectory")
# ------------------------------------------------------------------------------
