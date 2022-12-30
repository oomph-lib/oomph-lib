# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered picture_show subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(picture_show C CXX Fortran)

add_custom_target(
  generate_rotate_gifs_js ALL
  COMMAND ./turn_into_rotation_instructions.sh
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  BYPRODUCTS rotate_gifs.js)

message(VERBOSE "Leaving picture_show subdirectory")
# ------------------------------------------------------------------------------
