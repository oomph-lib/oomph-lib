# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered doc subdirectory")

# Project meta & global options
cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(doc LANGUAGES C CXX Fortran)

# This doc build is *configure + build* only
set(CMAKE_SKIP_INSTALL_RULES ON)

# Root of the whole repo
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/..")

# Generate the top-level doc from doc.txt
include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}" DOCFILE doc.txt)

# Demo-driver index pages
set(MAKE_INDEX_SCRIPT
    "${OOMPH_ROOT_DIR}/scripts/make_index_html_in_demo_driver_dirs.bash")

# Ask the helper script what by-products it will create
execute_process(
  COMMAND "${MAKE_INDEX_SCRIPT}" --list-byproducts
  OUTPUT_VARIABLE _index_raw
  OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE "\n" ";" INDEX_HTML_BYPRODUCTS "${_index_raw}")

add_custom_target(
  make_index_html_in_demo_drivers ALL
  COMMAND "${MAKE_INDEX_SCRIPT}"
  BYPRODUCTS ${INDEX_HTML_BYPRODUCTS}
  COMMENT "Generating index.html files for demo drivers")

# TODO: Update "the_distribution"; it's based on the soon-to-be-deprecated
# Makefile.am stuff

# TODO: Switch to using add_custom_command(...) in OomphGenerateDocFrom.cmake
# instead of add_custom_target(...) to avoid constantly rebuilding...

# Recurse!
set(SUBDIRS
    people
    index
    picture_show
    optimisation
    contact
    publications
    coding_conventions
    copyright
    intro
    FAQ
    meshes
    tar_files
    the_distribution
    example_code_list
    the_data_structure
    quick_guide
    poisson
    unsteady_heat
    linear_wave
    helmholtz
    navier_stokes
    axisym_navier_stokes
    advection_diffusion
    beam
    shell
    solid
    interaction
    multi_physics
    eigenproblems
    leftovers
    order_of_action_functions
    hijacking
    change_log
    bugs
    linear_solvers
    paraview
    preconditioners
    to_be_written
    young_laplace
    mpi
    creating_doc
    linear_elasticity
    time_harmonic_fourier_decomposed_linear_elasticity
    fourier_decomposed_helmholtz
    time_harmonic_linear_elasticity
    pml_fourier_decomposed_helmholtz
    pml_helmholtz
    axisym_linear_elasticity
    acoustic_fsi
    fourier_decomposed_acoustic_fsi
    pml_time_harmonic_linear_elasticity
    deliberately_broken_doc_for_checking
    acknowledgements
    search_results)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving doc subdirectory")
# ------------------------------------------------------------------------------
