# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_entry_flow subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(three_d_entry_flow C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE three_d_entry_flow.txt)

message(VERBOSE "Leaving three_d_entry_flow subdirectory")
# ------------------------------------------------------------------------------
