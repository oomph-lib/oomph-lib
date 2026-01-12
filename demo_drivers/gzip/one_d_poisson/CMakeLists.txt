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
message(VERBOSE "Entered one_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(one_d_poisson C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

set(EXTRA_LIBS "")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    list(APPEND EXTRA_LIBS stdc++fs)
  endif()
endif()

if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
  oomph_add_executable(
    NAME one_d_poisson_with_gzip_validation
    SOURCES one_d_poisson_with_gzip_validation.cc
    LIBRARIES oomph::poisson oomph::meshes oomph::generic ${EXTRA_LIBS})
endif()

# ---------------------------------- TESTING ----------------------------------

if(CMAKE_CXX_STANDARD MATCHES 17|20|23|26)
  oomph_define_test_data(validata validate.sh)
  oomph_add_pure_cpp_test(
    TEST_NAME gzip.one_d_poisson
    DEPENDS_ON one_d_poisson_with_gzip_validation
    LOG_FILE validation.log)
endif()

message(VERBOSE "Leaving one_d_poisson subdirectory")
# ------------------------------------------------------------------------------
