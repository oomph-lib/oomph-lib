# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_3d_ALE subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_adaptive_3d_ALE C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

if(OOMPH_HAS_CGAL)
  oomph_add_executable(
    NAME uns_adapt_3d
    SOURCES uns_adapt_3d.cc
    LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
    CXX_STANDARD 17)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_CGAL)
  oomph_add_test(
    TEST_NAME navier_stokes.unstructured_adaptive_3d_ALE
    TARGET_DEPENDENCIES uns_adapt_3d
    EXTRA_REQUIRES
    LABELS "navier_stokes;unstructured_adaptive_3d_ALE;adapt")
endif()

message(VERBOSE "Leaving unstructured_adaptive_3d_ALE subdirectory")
# ------------------------------------------------------------------------------
