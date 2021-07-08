#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
hp_adaptive_driven_cavity

#----------------------------------------------------------------------

# Sources for executable
hp_adaptive_driven_cavity_SOURCES = hp_adaptive_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hp_adaptive_driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


EXTRA_DIST += hp_adaptive_cavity_partition.dat


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

