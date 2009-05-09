#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
fsi_osc_ring_external

#----------------------------------------------------------------------

# Sources for executable
fsi_osc_ring_external_SOURCES = fsi_osc_ring.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_osc_ring_external_LDADD = \
                -L@libdir@ -lbeam -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
fsi_osc_ring_external_CXXFLAGS = -DUSE_EXTERNAL_ELEMENTS

#----------------------------------------------------------------------

EXTRA_DIST += fsi_osc_ring_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  

