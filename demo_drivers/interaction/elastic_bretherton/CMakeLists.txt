# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered elastic_bretherton subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(elastic_bretherton C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_bretherton
  SOURCES elastic_bretherton.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::fluid_interface
  CXX_STANDARD 14)

# Add a sanitiser string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}") # Create
# hash string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH) # Shorten to 7 characters
# set(HASHED_TARGET_NAME elastic_bretherton_${PATH_HASH}) # Append hash
# target_compile_options(${HASHED_TARGET_NAME} PRIVATE -fsanitize=address -g)
# target_link_options(${HASHED_TARGET_NAME} PRIVATE -fsanitize=address -g)
# if(NOT APPLE) target_compile_options(${HASHED_TARGET_NAME} PRIVATE
# -fsanitize=leak) target_link_options(${HASHED_TARGET_NAME} PRIVATE
# -fsanitize=leak) endif()

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.elastic_bretherton
  TARGET_DEPENDENCIES elastic_bretherton
  EXTRA_REQUIRES
  LABELS interaction elastic_bretherton)

message(VERBOSE "Leaving elastic_bretherton subdirectory")
# ------------------------------------------------------------------------------
