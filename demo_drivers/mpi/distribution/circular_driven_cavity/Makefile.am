#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
circular_driven_cavity

#----------------------------------------------------------------------

# Sources for executable
circular_driven_cavity_SOURCES = circular_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
circular_driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


EXTRA_DIST += circular_cavity_1_partition.dat \
              circular_cavity_2_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

