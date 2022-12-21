# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered orr_sommerfeld subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(orr_sommerfeld C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME orr_sommerfeld
    SOURCES orr_sommerfeld.cc
    LIBRARIES oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME eigenproblems.orr_sommerfeld
    TARGET_DEPENDENCIES orr_sommerfeld
    EXTRA_REQUIRES
    LABELS eigenproblems orr_sommerfeld)
endif()

message(VERBOSE "Leaving orr_sommerfeld subdirectory")
# ------------------------------------------------------------------------------
