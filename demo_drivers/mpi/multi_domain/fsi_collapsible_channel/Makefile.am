#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
fsi_collapsible_channel_adapt_external \
fsi_collapsible_channel_macro_adapt_external

#----------------------------------------------------------------------

# Sources for executable
fsi_collapsible_channel_adapt_external_SOURCES = fsi_collapsible_channel_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_adapt_external_LDADD = \
                -L@libdir@ -lbeam -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
fsi_collapsible_channel_adapt_external_CXXFLAGS = -DUSE_EXTERNAL_ELEMENTS

#----------------------------------------------------------------------

# Sources for executable
fsi_collapsible_channel_macro_adapt_external_SOURCES = fsi_collapsible_channel_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_collapsible_channel_macro_adapt_external_LDADD = \
                -L@libdir@ -lbeam -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
fsi_collapsible_channel_macro_adapt_external_CXXFLAGS = -DMACRO_ELEMENT_NODE_UPDATE -DUSE_EXTERNAL_ELEMENTS

#----------------------------------------------------------------------



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  

