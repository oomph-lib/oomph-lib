# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_ring subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unsteady_ring C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE unsteady_ring.txt)

message(VERBOSE "Leaving unsteady_ring subdirectory")
# ------------------------------------------------------------------------------
