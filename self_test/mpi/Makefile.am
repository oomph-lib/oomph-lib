# Do we want to run self tests which involves more than
# two processors?
if WANT_MPI_VARIABLENP
  variablenp_dir = \
		distribution_concatenation \
		vector_split_without_communication \
		vector_concatenation_without_communication \
		matrix_concatenation_without_communication \
		vector_split \
		vector_concatenation \
		matrix_concatenation \
		crdoublematrix_copy_constructor \
		matrix_addition 
else
  variablenp_dir = 
endif

SUBDIRS = \
$(variablenp_dir) \
two_d_mesh_dist \
three_d_mesh_dist \
line_visualiser




