# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_interface subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(adaptive_interface C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE adaptive_interface.txt)

message(VERBOSE "Leaving adaptive_interface subdirectory")
# ------------------------------------------------------------------------------
