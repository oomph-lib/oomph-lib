#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_fluid_spherical_cap 

# Sources for executable
two_fluid_spherical_cap_SOURCES = two_fluid_spherical_cap.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_fluid_spherical_cap_LDADD = -L@libdir@ -laxisym_navier_stokes \
				-lfluid_interface \
                                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
