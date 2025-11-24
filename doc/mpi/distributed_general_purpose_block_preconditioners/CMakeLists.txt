# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(
  VERBOSE
  "Entered distributed_general_purpose_block_preconditioners subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(distributed_general_purpose_block_preconditioners C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE distributed_general_purpose_block_preconditioners.txt)
message(
  VERBOSE
  "Leaving distributed_general_purpose_block_preconditioners subdirectory")
# ------------------------------------------------------------------------------
