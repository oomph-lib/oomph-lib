# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_dim subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_dim C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_dim
  SOURCES homo_lin_elasticity_elements.h homo_lin_elasticity_elements.cc
          Thomo_lin_elasticity_elements.h Thomo_lin_elasticity_elements.cc
          two_dim.cc
  LIBRARIES oomph::linear_elasticity
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_dim_hex
  SOURCES homo_lin_elasticity_elements.h homo_lin_elasticity_elements.cc
          Thomo_lin_elasticity_elements.h Thomo_lin_elasticity_elements.cc
          two_dim_hex.cc
  LIBRARIES oomph::linear_elasticity
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME homogenisation.two_dim
  TARGET_DEPENDENCIES two_dim two_dim_hex
  EXTRA_REQUIRES
  LABELS homogenisation two_dim)

message(VERBOSE "Leaving two_dim subdirectory")
# ------------------------------------------------------------------------------
