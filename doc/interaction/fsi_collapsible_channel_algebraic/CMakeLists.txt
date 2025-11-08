# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_collapsible_channel_algebraic subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_collapsible_channel_algebraic C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE fsi_chan_alg.txt)

message(VERBOSE "Leaving fsi_collapsible_channel_algebraic subdirectory")
# ------------------------------------------------------------------------------
