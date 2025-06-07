# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered generic subdirectory")

# Define the headers
set(HEADERS
    oomph_definitions.h
    Qelements.h
    Qspectral_elements.h
    elements.h
    Qelement_face_coordinate_translation_schemes.h
    integral.h
    assembly_handler.h
    periodic_orbit_handler.h
    problem.h
    linear_solver.h
    shape.h
    Vector.h
    frontal_solver.h
    matrices.h
    spines.h
    element_with_moving_nodes.h
    complex_matrices.h
    displacement_control_element.h
    mesh.h
    timesteppers.h
    explicit_timesteppers.h
    hermite_elements.h
    nodes.h
    oomph_utilities.h
    elastic_problems.h
    hijacked_elements.h
    geom_objects.h
    algebraic_elements.h
    macro_element.h
    stored_shape_function_elements.h
    map_matrix.h
    domain.h
    quad_mesh.h
    quadtree.h
    dg_elements.h
    error_estimator.h
    refineable_mesh.h
    fsi.h
    octree.h
    tree.h
    refineable_elements.h
    refineable_quad_element.h
    hp_refineable_elements.h
    refineable_line_spectral_element.h
    refineable_quad_spectral_element.h
    refineable_brick_spectral_element.h
    orthpoly.h
    refineable_quad_mesh.h
    pseudosolid_node_update_elements.h
    refineable_brick_element.h
    brick_mesh.h
    refineable_brick_mesh.h
    macro_element_node_update_element.h
    mesh_as_geometric_object.h
    face_element_as_geometric_object.h
    eigen_solver.h
    Telements.h
    triangle_scaffold_mesh.h
    geompack_scaffold_mesh.h
    tetgen_scaffold_mesh.h
    pseudo_buckling_ring.h
    simple_cubic_scaffold_tet_mesh.h
    line_mesh.h
    binary_tree.h
    refineable_line_element.h
    refineable_line_mesh.h
    triangle_mesh.h
    tet_mesh.h
    partitioning.h
    communicator.h
    linear_algebra_distribution.h
    double_vector.h
    double_multi_vector.h
    double_vector_with_halo.h
    multi_domain.h
    element_with_external_element.h
    iterative_linear_solver.h
    missing_masters.h
    preconditioner.h
    general_purpose_preconditioners.h
    block_preconditioner.h
    general_purpose_block_preconditioners.h
    SuperLU_preconditioner.h
    matrix_vector_product.h
    projection.h
    line_visualiser.h
    sum_of_matrices.h
    implicit_midpoint_rule.h
    trapezoid_rule.h
    preconditioner_array.h
    pml_meshes.h
    pml_mapping_functions.h
    generalised_timesteppers.h
    vector_matrix.h
    face_mesh_project.h
    generalised_newtonian_constitutive_models.h
    unstructured_two_d_mesh_geometry_base.h
    geometric_multigrid.h
    sample_point_container.h
    sample_point_parameters.h
    sparse_vector.h
    geom_obj_with_boundary.h
    extruded_macro_element.h
    extruded_domain.h
    black_box_newton_solver.h
    gzip_reader.h
    fpdiff.h)

if(OOMPH_HAS_STACKTRACE)
  list(APPEND HEADERS stacktrace.h)
endif()

# Define the sources
set(SOURCES
    oomph_definitions.cc
    oomph_utilities.cc
    complex_matrices.cc
    matrices.cc
    timesteppers.cc
    explicit_timesteppers.cc
    integral.cc
    nodes.cc
    elements.cc
    mesh.cc
    assembly_handler.cc
    periodic_orbit_handler.cc
    problem.cc
    Qelements.cc
    Qspectral_elements.cc
    frontal_solver.cc
    linear_solver.cc
    Qelement_face_coordinate_translation_schemes.cc
    Telements.cc
    hermite_elements.cc
    elastic_problems.cc
    hijacked_elements.cc
    algebraic_elements.cc
    macro_element.cc
    stored_shape_function_elements.cc
    quad_mesh.cc
    domain.cc
    quadtree.cc
    dg_elements.cc
    error_estimator.cc
    refineable_elements.cc
    pseudosolid_node_update_elements.cc
    refineable_quad_element.cc
    refineable_mesh.cc
    hp_refineable_elements.cc
    fsi.cc
    octree.cc
    tree.cc
    orthpoly.cc
    superlu.c
    superlu_complex.c
    refineable_brick_element.cc
    brick_mesh.cc
    spines.cc
    element_with_moving_nodes.cc
    macro_element_node_update_element.cc
    mesh_as_geometric_object.cc
    multi_domain.cc
    missing_masters.cc
    element_with_external_element.cc
    eigen_solver.cc
    triangle_scaffold_mesh.cc
    geompack_scaffold_mesh.cc
    tetgen_scaffold_mesh.cc
    simple_cubic_scaffold_tet_mesh.cc
    line_mesh.cc
    binary_tree.cc
    refineable_line_element.cc
    triangle_mesh.cc
    tet_mesh.cc
    partitioning.cc
    communicator.cc
    linear_algebra_distribution.cc
    double_vector.cc
    double_multi_vector.cc
    double_vector_with_halo.cc
    iterative_linear_solver.cc
    general_purpose_preconditioners.cc
    block_preconditioner.cc
    matrix_vector_product.cc
    sum_of_matrices.cc
    implicit_midpoint_rule.cc
    preconditioner_array.cc
    general_purpose_block_preconditioners.cc
    pml_meshes.cc
    unstructured_two_d_mesh_geometry_base.cc
    sample_point_container.cc
    sample_point_parameters.cc
    geometric_multigrid.cc
    extruded_macro_element.cc
    extruded_domain.cc
    black_box_newton_solver.cc
    gzip_reader.cc
    fpdiff.cc)

if(OOMPH_HAS_SUPERLU_DIST)
  list(APPEND SOURCES superlu_dist.c)
  set_source_files_properties(
    superlu_dist.c PROPERTIES COMPILE_FLAGS -Wno-implicit-function-declaration
                              -Wno-enum-conversion -Wno-switch)
endif()

if(OOMPH_HAS_HYPRE)
  list(APPEND HEADERS hypre_solver.h)
  list(APPEND SOURCES hypre_solver.cc)
endif()

if(OOMPH_HAS_TRILINOS)
  list(APPEND HEADERS trilinos_helpers.h trilinos_preconditioners.h
       trilinos_solver.h trilinos_eigen_solver.h)
  list(APPEND SOURCES trilinos_helpers.cc trilinos_preconditioners.cc
       trilinos_solver.cc trilinos_eigen_solver.cc)
endif()

if(OOMPH_HAS_MUMPS)
  list(APPEND HEADERS mumps_solver.h)
  list(APPEND SOURCES mumps_fortran_solver.F mumps_solver.cc)
endif()

# Need to worry about these for dependency but they shouldn't be included into
# the include directory
set(HEADERS_NO_INSTALL cfortran.h frontal.h lapack_qz.h mumps.h)

# Template-only files: Need to be installed to include directory but should not
# be listed in combined header file
set(TEMPLATES tree.template.cc refineable_mesh.template.cc
    multi_domain.template.cc missing_masters.template.cc)

# The generic library relies on the external sources
set(LINKLIBS ${EXTERNAL_SRC_LIBS})

# Set the name of the library
set(LIBNAME generic)

# ------------------------------------------------------------------------------
# NOTE: The files under TEMPLATES are specified as SOURCES_NO_BUILD to ensure
# they will be set as library dependencies but marked as HEADER_FILE_ONLY. As a
# result, the templates will be installed to build/include/ and will not be
# built into the generic library artefact. This is similar to what was done
# before with the Autotools framework but the TEMPLATES were not marked as
# library dependencies.
# ------------------------------------------------------------------------------
# Import the OomphLibraryConfig module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  SOURCES_NO_BUILD ${TEMPLATES}
  HEADERS_NO_INSTALL ${HEADERS_NO_INSTALL}
  LINKLIBS ${LINKLIBS})

# Headers in external_src are typically #include-ed by files in src/generic/
# like they sit in the generic folder, but that isn't the case at build-time. It
# is only true once the library has been installed because they'll all be
# bundled together in build/include/. To make things work at build time, we need
# to specify the location of the external_src headers.
target_include_directories(
  ${LIBNAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/external_src>)

# Silence some annoying warnings produced by this library
include(OomphTargetSilenceWarnings)
oomph_target_silence_warnings(
  TARGET ${LIBNAME}
  CXX_COMPILE_FLAGS
    -Wno-potentially-evaluated-expression
    -Wno-instantiation-after-specialization -Wno-tautological-compare)

# ------------------------------------------------------------------------------
# Handle third-party library stuff:
# ------------------------------------------------------------------------------
# If oomph-lib needs MPI, link it publically so that executables using the
# library also inherit MPI
if(OOMPH_HAS_MPI)
  target_link_libraries(${LIBNAME} PUBLIC MPI::MPI_C MPI::MPI_CXX
                                          MPI::MPI_Fortran)
endif()

target_link_libraries(${LIBNAME} PUBLIC ${OpenBLAS_LIBRARIES})

# NOTE: We publicly link to third-party libraries so that these libraries get
# propagated to the user. If a user pulls in the oomph-lib library they should
# automatically get access to the third-party libraries; they shouldn't have to
# manually link against all of those libraries themselves to use them.

# We need METIS for the partioning code. NOTE: There's subtlety to this choice
# of linking. We don't want to link generic directly to the METIS lib; we just
# need the metis.h header (in partitioning.h). We also want users who build
# source code/demo drivers with the generic library to be able to find metis.h
# when the generic library is used, so we have to make sure the linking type is
# PUBLIC. If we used PRIVATE linking, the sources of generic would inherit the
# include path to metis.h, but when we try to compile sources outside of the
# library, relying on generic, they wouldn't be able to find metis.h. If we use
# IMPORTED, the sources of generic won't inherit the include path. If we linked
# to METIS::METIS (the one that knows where the metis library is too) PUBLIC-ly,
# the linker would warn us about duplicate symbols because SuperLU provides the
# metis lib symbols. If we used INTERFACE linking instead, we wouldn't find the
# include directories again. That leaves PRIVATE linking, which we can't do with
# an "UNKNOWN IMPORTED" library.
target_link_libraries(${LIBNAME} PUBLIC METIS::METIS_HEADER_ONLY)

# Link to SuperLU. NOTE: SuperLU propagates its dependency on METIS
target_link_libraries(${LIBNAME} PUBLIC superlu::superlu)

if(OOMPH_HAS_SUPERLU_DIST)
  target_link_libraries(
    ${LIBNAME} PUBLIC SuperLU_DIST::SuperLU_DIST ParMETIS::ParMETIS
                      METIS::METIS GKlib::GKlib)
endif()

if(OOMPH_HAS_CGAL)
  target_link_libraries(${LIBNAME} PUBLIC CGAL::CGAL)
endif()

if(OOMPH_HAS_HYPRE)
  target_link_libraries(${LIBNAME} PUBLIC HYPRE::HYPRE)
endif()

if(OOMPH_HAS_TRILINOS)
  target_link_libraries(${LIBNAME} PUBLIC Trilinos::all_libs)
endif()

if(OOMPH_HAS_MUMPS)
  target_link_libraries(${LIBNAME} PUBLIC MUMPS::MUMPS)
endif()

message(VERBOSE "Leaving generic subdirectory")
# ------------------------------------------------------------------------------
