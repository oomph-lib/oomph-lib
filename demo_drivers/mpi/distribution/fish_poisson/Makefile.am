#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
fish_poisson

#----------------------------------------------------------------------

# Sources for executable
fish_poisson_SOURCES = fish_poisson.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------


EXTRA_DIST += fish_fully_automatic_partition.dat \
              fish_incremental_partition.dat \
              fish_selective_partition.dat distributed_adapt.mcr

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

