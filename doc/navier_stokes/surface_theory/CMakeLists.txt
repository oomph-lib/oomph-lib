# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered surface_theory subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(surface_theory C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE surface_theory.txt)

message(VERBOSE "Leaving surface_theory subdirectory")
# ------------------------------------------------------------------------------
