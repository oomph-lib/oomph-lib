# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered example_code_list subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(example_code_list C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE example_code_list.txt
  SUPPRESS_LATEX_IN_THIS_DIRECTORY)

message(VERBOSE "Leaving example_code_list subdirectory")
# ------------------------------------------------------------------------------
