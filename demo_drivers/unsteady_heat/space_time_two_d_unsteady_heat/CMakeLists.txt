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
message(VERBOSE "Entered space_time_two_d_unsteady_heat subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(space_time_two_d_unsteady_heat C CXX Fortran)

if(NOT oomphlib_FOUND)
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
  DEPENDS_ON test_mixed_order_galerkin_petrov test_equal_order_galerkin_petrov
             test_equal_order_galerkin
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving space_time_two_d_unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
