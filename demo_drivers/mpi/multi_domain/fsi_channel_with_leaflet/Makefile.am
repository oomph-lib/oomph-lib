#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
fsi_channel_with_leaflet

#----------------------------------------------------------------------

# Sources for executable
fsi_channel_with_leaflet_SOURCES = fsi_channel_with_leaflet.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_channel_with_leaflet_LDADD = \
                -L@libdir@ -lbeam -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

EXTRA_DIST += fsi_channel_with_leaflet_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

