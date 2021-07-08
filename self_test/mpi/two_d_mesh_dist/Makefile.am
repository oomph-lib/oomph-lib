#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
two_d_mesh_dist

#----------------------------------------------------------------------

# Sources for executable
two_d_mesh_dist_SOURCES = two_d_mesh_dist.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_mesh_dist_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

EXTRA_DIST += two_d_mesh_dist_partition.dat

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

