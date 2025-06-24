# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_3d_ALE subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_adaptive_3d_ALE C CXX Fortran)

if(NOT oomphlib_FOUND)
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
    DEPENDS_ON uns_adapt_3d
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
    TEST_FILES validate.sh validata)
endif()

message(VERBOSE "Leaving unstructured_adaptive_3d_ALE subdirectory")
# ------------------------------------------------------------------------------
