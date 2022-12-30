# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered doc subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(doc C CXX Fortran)

set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../")
include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}" DOCFILE doc.txt)

# Target to build documentation for all subdirectories of doc; dependencies will
# be added to this target through calls to oomph_generate_doc_from(...) in the
# subdirectories
add_custom_target(build_docs ALL)

# Get the list of index.html files that will be produced by
# make_index_html_in_demo_driver_dirs.bash
execute_process(
  COMMAND "${OOMPH_ROOT_DIR}/scripts/make_index_html_in_demo_driver_dirs.bash"
          --list-byproducts OUTPUT_VARIABLE INDEX_HTML_BYPRODUCTS)
string(REPLACE "\n" ";" INDEX_HTML_BYPRODUCTS_AS_LIST ${INDEX_HTML_BYPRODUCTS})

# Target to generate index.html files for the demo drivers
add_custom_target(
  make_index_html_in_demo_drivers ALL
  COMMAND "${OOMPH_ROOT_DIR}/scripts/make_index_html_in_demo_driver_dirs.bash"
  BYPRODUCTS ${INDEX_HTML_BYPRODUCTS_AS_LIST})

set(SUBDIRS people)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving doc subdirectory")
# ------------------------------------------------------------------------------
