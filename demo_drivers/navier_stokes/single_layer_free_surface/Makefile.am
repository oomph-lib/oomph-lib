#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=single_layer 

# Sources for executable
single_layer_SOURCES = single_layer.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
single_layer_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)
