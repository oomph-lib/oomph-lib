#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
soluble_surfactant 


#---------------------------------------------------------------------


# Sources for executable
soluble_surfactant_SOURCES = double_buoyant_navier_stokes_elements.h soluble_surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
soluble_surfactant_LDADD = -L@libdir@ -ladvection_diffusion_reaction \
 -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
