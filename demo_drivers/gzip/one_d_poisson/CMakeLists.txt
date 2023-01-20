# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(one_d_poisson C CXX Fortran)
message("CMAKE_CXX_STANDARD: ${CMAKE_CXX_STANDARD}")
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

# FIXME: Compiler (on macOS) producing a warning about
# std::experimental::filesystem being deprecated into favour of C++17's
# std::filesystem and "This experimental version will be removed in LLVM 11.".
# Should revert to C++17 but will need to patch the issue with Trilinos if it is
# still around
if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
  oomph_add_executable(
    NAME one_d_poisson_with_gzip_validation
    SOURCES one_d_poisson_with_gzip_validation.cc
    LIBRARIES oomph::poisson oomph::meshes oomph::generic stdc++fs)
endif()

# ---------------------------------- TESTING ----------------------------------

if(CMAKE_CXX_STANDARD MATCHES 17|20|23|26)
  oomph_define_test_data(validata validate.sh)
  oomph_add_pure_cpp_test(
    TEST_NAME gzip.one_d_poisson
    DEPENDS_ON one_d_poisson_with_gzip_validation
    LOG_FILE validation.log
    LABELS gzip one_d_poisson)
endif()

message(VERBOSE "Leaving one_d_poisson subdirectory")
# ------------------------------------------------------------------------------
