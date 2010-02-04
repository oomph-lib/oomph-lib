#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_layer_interface_axisym

# Sources for executable
two_layer_interface_axisym_SOURCES = two_layer_interface_axisym.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
two_layer_interface_axisym_LDADD = -L@libdir@ -laxisym_navier_stokes -lnavier_stokes -lfluid_interface -lgeneric -loomph_crbond_bessel $(EXTERNAL_LIBS) $(FLIBS)
