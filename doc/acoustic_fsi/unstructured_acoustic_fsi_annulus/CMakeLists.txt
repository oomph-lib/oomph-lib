# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_acoustic_fsi_annulus subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(unstructured_acoustic_fsi_annulus C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE unstructured_acoustic_fsi_annulus.txt)

message(VERBOSE "Leaving unstructured_acoustic_fsi_annulus subdirectory")
# ------------------------------------------------------------------------------
