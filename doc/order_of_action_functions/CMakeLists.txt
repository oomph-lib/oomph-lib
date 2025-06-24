# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered order_of_action_functions subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(order_of_action_functions C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE order_of_action_functions.txt
  SUPPRESS_LATEX_IN_THIS_DIRECTORY)

message(VERBOSE "Leaving order_of_action_functions subdirectory")
# ------------------------------------------------------------------------------
