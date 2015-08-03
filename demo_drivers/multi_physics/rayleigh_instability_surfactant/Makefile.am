#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
rayleigh_instability_insoluble_surfactant 


#---------------------------------------------------------------------


# Sources for executable
rayleigh_instability_insoluble_surfactant_SOURCES = rayleigh_instability_insoluble_surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rayleigh_instability_insoluble_surfactant_LDADD = -L@libdir@ \
 -laxisym_navier_stokes -lmulti_physics -lfluid_interface -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)
