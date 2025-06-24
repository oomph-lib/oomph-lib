# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered contact_angle subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(contact_angle C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE contact_angle.txt)

message(VERBOSE "Leaving contact_angle subdirectory")
# ------------------------------------------------------------------------------
