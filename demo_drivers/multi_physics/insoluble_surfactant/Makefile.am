#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
surfactant 


#---------------------------------------------------------------------


# Sources for executable
surfactant_SOURCES = surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
surfactant_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lmulti_physics -lfluid_interface -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)
