#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=spherical_shell_convection

# Sources for executable
spherical_shell_convection_SOURCES = spherical_shell_convection.cc \
                                 spherical_buoyant_navier_stokes.h
# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spherical_shell_convection_LDADD = -L@libdir@ \
 -lspherical_advection_diffusion \
 -lspherical_navier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
