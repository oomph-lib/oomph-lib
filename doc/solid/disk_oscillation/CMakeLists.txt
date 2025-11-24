# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disk_oscillation subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(disk_oscillation C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE disk_oscillation.txt)

message(VERBOSE "Leaving disk_oscillation subdirectory")
# ------------------------------------------------------------------------------
