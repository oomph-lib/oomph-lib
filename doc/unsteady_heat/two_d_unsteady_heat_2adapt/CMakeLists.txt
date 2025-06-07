# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_2adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_unsteady_heat_2adapt C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE two_d_unsteady_heat_2adapt.txt)

# The symlink we want to create
set(SRC_GIF
    "${CMAKE_CURRENT_SOURCE_DIR}/../two_d_unsteady_heat_ALE/figures/domain.gif")
set(DST_GIF "${CMAKE_CURRENT_SOURCE_DIR}/figures/domain.gif")

# Make sure we have a symlink
add_custom_command(
  OUTPUT "${DST_GIF}"
  COMMAND
    ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ensure_symlink.cmake"
    "${SRC_GIF}" "${DST_GIF}"
  DEPENDS "${SRC_GIF}"
  COMMENT "Ensuring symlink ${DST_GIF} --> ${SRC_GIF}"
  VERBATIM)

# Define a target to run the above command in the build
add_custom_target(generate_domain_gif ALL DEPENDS "${DST_GIF}")

# Don't delete the symlink when the 'clean' target is run
set_source_files_properties("${DST_GIF}" PROPERTIES SKIP_CLEAN TRUE)

message(VERBOSE "Leaving two_d_unsteady_heat_2adapt subdirectory")
# ------------------------------------------------------------------------------
