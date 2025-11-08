# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered deliberately_broken_doc_for_checking subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(deliberately_broken_doc_for_checking C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE deliberately_broken_doc_for_checking.txt)

message(VERBOSE "Leaving deliberately_broken_doc_for_checking subdirectory")
# ------------------------------------------------------------------------------
