#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
adaptive_driven_cavity \
adaptive_driven_cavity_load_balance

#----------------------------------------------------------------------

# Sources for executable
adaptive_driven_cavity_SOURCES = adaptive_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
adaptive_driven_cavity_load_balance_SOURCES = \
 adaptive_driven_cavity_load_balance.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_driven_cavity_load_balance_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


EXTRA_DIST += adaptive_cavity_1_partition.dat \
              adaptive_cavity_2_partition.dat \
              load_balance.mcr prune_and_balance.mcr


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

