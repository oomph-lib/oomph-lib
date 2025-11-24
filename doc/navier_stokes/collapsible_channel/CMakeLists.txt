# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(collapsible_channel C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE collapsible_channel.txt)

message(VERBOSE "Leaving collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
