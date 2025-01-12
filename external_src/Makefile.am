# Subdirectories get built into individual libraries
# Note that we specify . explicitly otherwise we descend into
# subdirectories first; see
#
#   https://www.gnu.org/software/automake/manual/html_node/Subdirectories.html
#
SUBDIRS = . oomph_hsl oomph_crbond_bessel 


# Build malloc.h in parmetis directory if needed (do it from here and
# before descending into subdir because it may be build after
# library itself (in which case the compilation will obviously fail.
if ALREADY_HAVE_METIS
else
all-am:
	cd oomph_parmetis_3.1.1; make dummy; pwd; ls -l malloc*; cd ..
endif

# Suppress triangle?
if SUPPRESS_TRIANGLE_LIB
else
  SUBDIRS += oomph_triangle
endif

# Suppress tetgen?
if SUPPRESS_TETGEN_LIB
else
  SUBDIRS += oomph_tetgen
endif


# Check if metis is already available
if ALREADY_HAVE_METIS
else
#  SUBDIRS += oomph_metis_4.0 
  SUBDIRS += oomph_metis_from_parmetis_3.1.1
endif

# Check if superlu is already available
if ALREADY_HAVE_SUPERLU
else
  SUBDIRS += oomph_superlu_4.3 
endif


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

# Include SuperLU_dist and parmetis if compiled with MPI
if WANT_MPI

# Check if superlu is already available
if ALREADY_HAVE_SUPERLU_DIST
else
  SUBDIRS += oomph_superlu_dist_3.0 
endif


# Check if parmetis is already available
if ALREADY_HAVE_PARMETIS
else
  SUBDIRS += oomph_parmetis_3.1.1 
endif

endif

# Tidy up: Remove symbolic links that make the combined header files appear
# in the directory above the actual headers as in the installed include
# directory.
clean-local:
	rm -f oomph_metis_4.0.h oomph_superlu_4.3.h oomph_parmetis_3.1.1.h


EXTRA_DIST=parmetis.h
