# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered raw_doxygen subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(raw_doxygen C CXX Fortran)

add_custom_target(
  build_raw_doxygen ALL
  COMMAND doxygen
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Files to clean up with "ninja clean" or "make clean"; can't just specify html/
# via BYPRODUCTS as its a directory and Ninja/Makefile generators won't wipe
# them if they're populated. They will if you use set_property(...) though
set_property(
  TARGET build_raw_doxygen
  APPEND
  PROPERTY ADDITIONAL_CLEAN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/html/")

message(VERBOSE "Leaving raw_doxygen subdirectory")
# ------------------------------------------------------------------------------
