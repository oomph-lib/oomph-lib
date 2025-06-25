# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered lgr_navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(lgr_navier_stokes C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE lgr_navier_stokes.txt)

message(VERBOSE "Leaving lgr_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
