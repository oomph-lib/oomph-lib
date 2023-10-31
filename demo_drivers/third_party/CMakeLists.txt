# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered third_party subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(third_party C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

set(SUBDIRS)

if(OOMPH_HAS_NLOHMANN_JSON)
  list(APPEND SUBDIRS test_nlohmann_json)
endif()

# Add desired subdirectories to the build
foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving third_party subdirectory")
# ------------------------------------------------------------------------------
