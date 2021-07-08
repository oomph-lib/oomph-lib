#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS=
#Only check if we have Trilinos
if OOMPH_HAS_TRILINOS
check_PROGRAMS += harmonic
endif

#----------------------------------------------------------------------

# Sources for executable
harmonic_SOURCES = harmonic.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
harmonic_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


#EXTRA_DIST += harmonic_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

