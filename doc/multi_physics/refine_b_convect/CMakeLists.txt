# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered refine_b_convect subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(refine_b_convect C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE refine_b_convect.txt)

message(VERBOSE "Leaving refine_b_convect subdirectory")
# ------------------------------------------------------------------------------
