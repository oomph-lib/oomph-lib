# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered adaptive_driven_cavity_load_balance subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(adaptive_driven_cavity_load_balance C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE adaptive_driven_cavity_load_balance.txt)
message(VERBOSE "Leaving adaptive_driven_cavity_load_balance subdirectory")
# ------------------------------------------------------------------------------
