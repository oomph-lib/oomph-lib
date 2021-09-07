
# A few file definitions
#-----------------------



# Define the sources
sources =  \
oomph_definitions.cc oomph_utilities.cc \
complex_matrices.cc \
matrices.cc       timesteppers.cc explicit_timesteppers.cc \
integral.cc   nodes.cc  \
elements.cc mesh.cc assembly_handler.cc periodic_orbit_handler.cc problem.cc \
Qelements.cc Qspectral_elements.cc frontal_solver.cc      linear_solver.cc \
Qelement_face_coordinate_translation_schemes.cc \
Telements.cc hermite_elements.cc   \
elastic_problems.cc  hijacked_elements.cc      \
algebraic_elements.cc \
macro_element.cc \
stored_shape_function_elements.cc \
quad_mesh.cc \
domain.cc   quadtree.cc \
dg_elements.cc \
error_estimator.cc   \
refineable_elements.cc pseudosolid_node_update_elements.cc \
refineable_quad_element.cc  refineable_mesh.cc \
hp_refineable_elements.cc \
fsi.cc octree.cc tree.cc orthpoly.cc             \
superlu.c superlu_complex.c refineable_brick_element.cc \
brick_mesh.cc spines.cc element_with_moving_nodes.cc \
macro_element_node_update_element.cc \
mesh_as_geometric_object.cc \
multi_domain.cc \
missing_masters.cc \
element_with_external_element.cc \
eigen_solver.cc \
triangle_scaffold_mesh.cc  geompack_scaffold_mesh.cc \
tetgen_scaffold_mesh.cc simple_cubic_scaffold_tet_mesh.cc \
line_mesh.cc binary_tree.cc refineable_line_element.cc \
triangle_mesh.cc tet_mesh.cc \
partitioning.cc communicator.cc linear_algebra_distribution.cc \
double_vector.cc \
double_multi_vector.cc \
double_vector_with_halo.cc \
iterative_linear_solver.cc \
general_purpose_preconditioners.cc block_preconditioner.cc \
matrix_vector_product.cc \
sum_of_matrices.cc \
implicit_midpoint_rule.cc \
preconditioner_array.cc general_purpose_block_preconditioners.cc pml_meshes.cc \
unstructured_two_d_mesh_geometry_base.cc sample_point_container.cc \
sample_point_parameters.cc geometric_multigrid.cc \
extruded_macro_element.cc extruded_domain.cc \
black_box_newton_solver.cc

if OOMPH_HAS_MUMPS
sources+=mumps_fortran_solver.F mumps_solver.cc
endif

if WANT_MPI
sources += superlu_dist.c
endif


# Define the headers
headers =  \
oomph_definitions.h Qelements.h    Qspectral_elements.h        elements.h    \
Qelement_face_coordinate_translation_schemes.h \
integral.h       assembly_handler.h periodic_orbit_handler.h problem.h \
linear_solver.h       shape.h \
Vector.h            frontal_solver.h      matrices.h       spines.h \
element_with_moving_nodes.h \
complex_matrices.h \
displacement_control_element.h \
mesh.h           timesteppers.h  explicit_timesteppers.h \
hermite_elements.h  nodes.h      oomph_utilities.h \
elastic_problems.h  hijacked_elements.h      geom_objects.h \
algebraic_elements.h            macro_element.h \
stored_shape_function_elements.h \
map_matrix.h \
domain.h                        quad_mesh.h \
quadtree.h \
dg_elements.h \
error_estimator.h \
refineable_mesh.h \
fsi.h octree.h  tree.h  \
refineable_elements.h refineable_quad_element.h \
hp_refineable_elements.h \
refineable_line_spectral_element.h \
refineable_quad_spectral_element.h \
refineable_brick_spectral_element.h orthpoly.h \
refineable_quad_mesh.h pseudosolid_node_update_elements.h \
refineable_brick_element.h brick_mesh.h refineable_brick_mesh.h  \
macro_element_node_update_element.h \
mesh_as_geometric_object.h face_element_as_geometric_object.h \
eigen_solver.h \
Telements.h \
triangle_scaffold_mesh.h geompack_scaffold_mesh.h tetgen_scaffold_mesh.h \
pseudo_buckling_ring.h simple_cubic_scaffold_tet_mesh.h \
line_mesh.h binary_tree.h refineable_line_element.h \
refineable_line_mesh.h \
triangle_mesh.h tet_mesh.h \
partitioning.h communicator.h linear_algebra_distribution.h double_vector.h \
double_multi_vector.h double_vector_with_halo.h \
multi_domain.h element_with_external_element.h iterative_linear_solver.h \
missing_masters.h \
preconditioner.h \
general_purpose_preconditioners.h block_preconditioner.h \
general_purpose_block_preconditioners.h SuperLU_preconditioner.h \
matrix_vector_product.h projection.h line_visualiser.h \
sum_of_matrices.h implicit_midpoint_rule.h \
trapezoid_rule.h \
preconditioner_array.h pml_meshes.h pml_mapping_functions.h \
generalised_timesteppers.h vector_matrix.h face_mesh_project.h \
generalised_newtonian_constitutive_models.h \
unstructured_two_d_mesh_geometry_base.h \
geometric_multigrid.h sample_point_container.h \
sample_point_parameters.h sparse_vector.h \
geom_obj_with_boundary.h extruded_macro_element.h extruded_domain.h \
black_box_newton_solver.h


if OOMPH_HAS_MUMPS
headers+=mumps_solver.h
endif

if HAVE_STACKTRACE
headers += stacktrace.h
endif

# Template-only files: Need to be included in include directory
# but should not be listed in combined header file
templates = \
tree.template.cc refineable_mesh.template.cc multi_domain.template.cc \
missing_masters.template.cc

# Do we have the HYPRE library?
if OOMPH_HAS_HYPRE

   # Add the hypre sources
   sources += hypre_solver.cc


   # Add the hypre headers
   headers += hypre_solver.h

endif

# Do we have the TRILINOS library?
if OOMPH_HAS_TRILINOS

   # Add the trilinos sources
   sources += trilinos_helpers.cc trilinos_preconditioners.cc \
              trilinos_solver.cc trilinos_eigen_solver.cc


   # Add the trilinos headers
   headers +=   trilinos_helpers.h trilinos_preconditioners.h \
                trilinos_solver.h trilinos_eigen_solver.h

endif

# Need to worry about these for dependency but they shouldn't
# be included into the include directory
extra_headers = cfortran.h frontal.h arpack.h lapack_qz.h mumps.h

# Define name of library
libname = generic

# Combine headers and sources
headers_and_sources = $(headers) $(extra_headers) $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libgeneric.la


# Sources that the library depends on:
#-------------------------------------
libgeneric_la_SOURCES = $(headers_and_sources)


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/generic

#These are the header files that are to be placed in subdirectory
library_include_HEADERS=$(headers) $(templates)

# Required libraries
#-------------------
# Include path for library headers -- need to refer to include files
# in src as they will not have been installed yet!
AM_CPPFLAGS += -I$(top_srcdir)/external_src -I$(top_srcdir)/external_src/oomph_tetgen



# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
$(libname).h: dummy_$(libname).h

dummy_$(libname).h:  $(headers)
	echo  $(libname) $(headers) > all_$(libname).aux
	$(AWK) -f $(top_srcdir)/bin/headers.awk < \
	       all_$(libname).aux > $(libname).h
	rm all_$(libname).aux


# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook:
	(cd $(library_includedir) && rm -f $(headers) $(templates) )
	(echo "$(headers) $(templates)" > include_files.list )
	($(top_srcdir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list )
else
install-data-hook:
endif


# Cleanup: Kill automatically generated combined header file
clean-local:
	rm -f $(libname).h *~
distclean-local:
	rm -f $(libname).h *~

#EXTRA_DIST=prettyprint.h prettyprint98.h
