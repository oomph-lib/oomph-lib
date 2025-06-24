# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered index subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(index C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../")

# Get relative path to oomph-lib root
cmake_path(
  RELATIVE_PATH OOMPH_ROOT_DIR BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  OUTPUT_VARIABLE RELATIVE_PATH_TO_OOMPH_ROOT_DIR)

# Define a rule to generate oomph_index.txt; oomph_generate_doc_from(...) will
# generate a target that will make sure the docfile is produced before building
# the docs
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/oomph_index.txt"
  COMMAND cp intro.txt oomph_index.txt
  COMMAND "${OOMPH_ROOT_DIR}/scripts/make_index.py" oomph.idx
          "${RELATIVE_PATH_TO_OOMPH_ROOT_DIR}" >> oomph_index.txt
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Define target to build docs
include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE oomph_index.txt
  BUILD_DOXY_HEADER_TARGET)

# Strip trailing slash (not necessary, just for prettiness...)
string(REGEX REPLACE "/$" "" RELATIVE_PATH_TO_OOMPH_ROOT_DIR
                     "${RELATIVE_PATH_TO_OOMPH_ROOT_DIR}")

message(VERBOSE "Leaving index subdirectory")
# ------------------------------------------------------------------------------
