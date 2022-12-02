# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered space_time_two_d_unsteady_heat subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(space_time_two_d_unsteady_heat C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME test_mixed_order_galerkin_petrov
  SOURCES test_mixed_order_galerkin_petrov.cc
  LIBRARIES oomph::space_time_unsteady_heat_mixed_order_galerkin_petrov
            oomph::space_time_block_preconditioner oomph::meshes oomph::generic)

oomph_add_executable(
  NAME test_equal_order_galerkin_petrov
  SOURCES test_equal_order_galerkin_petrov.cc
  LIBRARIES oomph::space_time_unsteady_heat_equal_order_galerkin_petrov
            oomph::space_time_block_preconditioner oomph::meshes oomph::generic)

oomph_add_executable(
  NAME test_equal_order_galerkin
  SOURCES test_equal_order_galerkin.cc
  LIBRARIES oomph::space_time_unsteady_heat_equal_order_galerkin
            oomph::space_time_block_preconditioner oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME unsteady_heat.space_time_two_d_unsteady_heat
  TARGET_DEPENDENCIES test_mixed_order_galerkin_petrov
                      test_equal_order_galerkin_petrov test_equal_order_galerkin
  EXTRA_REQUIRES
  LABELS unsteady_heat space_time_two_d_unsteady_heat)

message(VERBOSE "Leaving space_time_two_d_unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
