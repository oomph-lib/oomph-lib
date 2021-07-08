#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=single_layer_free_surface_axisym

# Sources for executable
single_layer_free_surface_axisym_SOURCES = single_layer_free_surface_axisym.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
single_layer_free_surface_axisym_LDADD = -L@libdir@ -laxisym_navier_stokes -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
