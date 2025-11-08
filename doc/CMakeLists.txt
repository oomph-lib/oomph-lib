# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered doc subdirectory")

# Project meta & global options
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(doc LANGUAGES C CXX Fortran)

# This doc build is *configure + build* only
set(CMAKE_SKIP_INSTALL_RULES ON)

# Root of the whole repo
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/..")

# Generate the top-level doc from doc.txt
include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}" DOCFILE doc.txt)

# ------------------------------------------------------------------------------
# Demo-driver index pages
set(MAKE_INDEX_SCRIPT
    "${OOMPH_ROOT_DIR}/scripts/make_index_html_in_demo_driver_dirs.bash")

# Ask the helper script which files it will (re)create
execute_process(
  COMMAND "${MAKE_INDEX_SCRIPT}" --list-byproducts
  OUTPUT_VARIABLE _index_raw
  OUTPUT_STRIP_TRAILING_WHITESPACE)

# Convert its newline output into a true CMake list.
string(REPLACE "\n" ";" INDEX_HTML_FILES "${_index_raw}")

# One custom command that claims *every* html file as an OUTPUT
add_custom_command(
  OUTPUT ${INDEX_HTML_FILES}
  COMMAND "${MAKE_INDEX_SCRIPT}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  COMMENT "Generating index.html files for demo drivers"
  VERBATIM)

# A convenience target that the build can depend on
add_custom_target(make_index_html_in_demo_drivers ALL
                  DEPENDS ${INDEX_HTML_FILES})
# ------------------------------------------------------------------------------

# Recurse!
set(SUBDIRS
    FAQ
    acknowledgements
    acoustic_fsi
    advection_diffusion
    axisym_linear_elasticity
    axisym_navier_stokes
    beam
    bugs
    change_log
    coding_conventions
    contact
    copyright
    creating_doc
    deliberately_broken_doc_for_checking
    eigenproblems
    example_code_list
    fourier_decomposed_acoustic_fsi
    fourier_decomposed_helmholtz
    helmholtz
    hijacking
    index
    interaction
    intro
    leftovers
    linear_elasticity
    linear_solvers
    linear_wave
    meshes
    mpi
    multi_physics
    navier_stokes
    optimisation
    order_of_action_functions
    paraview
    people
    picture_show
    pml_fourier_decomposed_helmholtz
    pml_helmholtz
    pml_time_harmonic_linear_elasticity
    poisson
    preconditioners
    publications
    quick_guide
    search_results
    shell
    solid
    tar_files
    the_data_structure
    time_harmonic_fourier_decomposed_linear_elasticity
    time_harmonic_linear_elasticity
    to_be_written
    unsteady_heat
    young_laplace)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving doc subdirectory")
# ------------------------------------------------------------------------------
