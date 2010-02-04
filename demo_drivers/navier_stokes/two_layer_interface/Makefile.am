#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_layer_interface 

# Sources for executable
two_layer_interface_SOURCES = two_layer_interface.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
two_layer_interface_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
