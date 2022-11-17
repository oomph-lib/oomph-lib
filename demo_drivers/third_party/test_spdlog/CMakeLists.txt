# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered test_spdlog subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(test_spdlog C CXX Fortran)

if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME test_spdlog
  SOURCES test_spdlog.cc
  LIBRARIES spdlog::spdlog
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME third_party.test_spdlog
  TARGET_DEPENDENCIES test_spdlog
  EXTRA_REQUIRES
  LABELS third_party test_spdlog
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving test_spdlog subdirectory")
# ------------------------------------------------------------------------------
