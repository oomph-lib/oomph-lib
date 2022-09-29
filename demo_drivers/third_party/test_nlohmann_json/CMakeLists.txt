# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered test_nlohmann_json subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(test_nlohmann_json C CXX Fortran)

if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME test_nlohmann_json
  SOURCES test_nlohmann_json.cc
  LIBRARIES nlohmann_json::nlohmann_json
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME third_party.test_nlohmann_json
  TARGET_DEPENDENCIES test_nlohmann_json
  EXTRA_REQUIRES
  LABELS "third_party;test_nlohmann_json")

message(VERBOSE "Leaving test_nlohmann_json subdirectory")
# ------------------------------------------------------------------------------
