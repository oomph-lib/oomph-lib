# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered compressed_square subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(compressed_square C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE compressed_square.txt)

message(VERBOSE "Leaving compressed_square subdirectory")
# ------------------------------------------------------------------------------
