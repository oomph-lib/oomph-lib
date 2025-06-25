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
message(VERBOSE "Entered optimisation subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(optimisation C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

set(SUBDIRS
    C_style_output
    disable_ALE
    linear_vs_nonlinear
    stored_shape_fcts
    use_coarse_base_meshes
    sparse_assemble
    fsi_jacobian_approximation)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving optimisation subdirectory")
# ------------------------------------------------------------------------------
