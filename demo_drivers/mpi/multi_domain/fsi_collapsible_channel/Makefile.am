#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
fsi_collapsible_channel_adapt \
fsi_collapsible_channel_macro_adapt

#----------------------------------------------------------------------

# Sources for executable
fsi_collapsible_channel_adapt_SOURCES = fsi_collapsible_channel_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_adapt_LDADD = \
                -L@libdir@ -lbeam -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
fsi_collapsible_channel_macro_adapt_SOURCES = fsi_collapsible_channel_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_macro_adapt_LDADD = \
                -L@libdir@ -lbeam -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
fsi_collapsible_channel_macro_adapt_CXXFLAGS = -DMACRO_ELEMENT_NODE_UPDATE

#----------------------------------------------------------------------

EXTRA_DIST += fsi_collapsible_channel_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

