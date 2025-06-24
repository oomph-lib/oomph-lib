# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered osc_ellipse subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(osc_ellipse C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE osc_ellipse.txt)

message(VERBOSE "Leaving osc_ellipse subdirectory")
# ------------------------------------------------------------------------------
