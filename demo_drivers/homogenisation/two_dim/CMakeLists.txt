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
message(VERBOSE "Entered two_dim subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_dim C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME two_dim
  SOURCES homo_lin_elasticity_elements.h homo_lin_elasticity_elements.cc
          Thomo_lin_elasticity_elements.h Thomo_lin_elasticity_elements.cc
          two_dim.cc
  LIBRARIES oomph::linear_elasticity oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_dim_hex
  SOURCES homo_lin_elasticity_elements.h homo_lin_elasticity_elements.cc
          Thomo_lin_elasticity_elements.h Thomo_lin_elasticity_elements.cc
          two_dim_hex.cc
  LIBRARIES oomph::linear_elasticity oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME homogenisation.two_dim
  DEPENDS_ON two_dim two_dim_hex
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving two_dim subdirectory")
# ------------------------------------------------------------------------------
