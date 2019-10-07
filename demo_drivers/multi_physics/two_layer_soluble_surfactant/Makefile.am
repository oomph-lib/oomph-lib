#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
two_layer_soluble_surfactant #refineable_two_layer_soluble_surfactant


#---------------------------------------------------------------------


# Sources for executable
two_layer_soluble_surfactant_SOURCES = double_buoyant_navier_stokes_elements.h \
two_layer_soluble_surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_layer_soluble_surfactant_LDADD = -L@libdir@ -ladvection_diffusion_reaction \
 -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
#refineable_two_layer_soluble_surfactant_SOURCES = double_buoyant_navier_stokes_elements.h \
#refineable_two_layer_soluble_surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
#refineable_two_layer_soluble_surfactant_LDADD = -L@libdir@ -ladvection_diffusion_reaction \
# -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
