#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
rayleigh_instability_insoluble_surfactant \
3d_rayleigh_instability_surfactant \
rayleigh_instability_soluble_surfactant


#---------------------------------------------------------------------


# Sources for executable
rayleigh_instability_insoluble_surfactant_SOURCES = rayleigh_instability_insoluble_surfactant.cc axisymmetric_advection_navier_stokes_elements.h


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rayleigh_instability_insoluble_surfactant_LDADD = -L@libdir@ \
 -laxisym_navier_stokes -lfluid_interface -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
rayleigh_instability_soluble_surfactant_SOURCES = rayleigh_instability_soluble_surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rayleigh_instability_soluble_surfactant_LDADD = -L@libdir@ \
 -laxisym_navier_stokes -lfluid_interface -lsolid -lconstitutive -laxisym_advection_diffusion  -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
3d_rayleigh_instability_surfactant_SOURCES = 3d_rayleigh_instability_surfactant.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
3d_rayleigh_instability_surfactant_LDADD = -L@libdir@ \
 -lnavier_stokes -lfluid_interface -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)
