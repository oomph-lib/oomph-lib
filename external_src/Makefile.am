# Subdirectories get built into individual libraries
SUBDIRS = oomph_superlu_3.0 oomph_metis_4.0 oomph_hsl oomph_crbond_bessel oomph_arpack


# Check if blas is already available
if ALREADY_HAVE_BLAS
else
  SUBDIRS += oomph_blas 
endif

# Check if lapack is already available
if ALREADY_HAVE_LAPACK
else
  SUBDIRS += oomph_lapack oomph_flapack
endif

# Include SuperLU_dist if compiled with MPI
if WANT_MPI
  SUBDIRS += oomph_superlu_dist_2.0
endif

# Tidy up: Remove symbolic links that make the combined header files appear
# in the directory above the actual headers as in the installed include
# directory.
clean-local:
	rm -f oomph_metis_4.0.h oomph_superlu_3.0.h oomph_superlu_dist_2.0.h
