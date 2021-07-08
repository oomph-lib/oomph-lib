#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=rayleigh_channel rayleigh_circular_driven_cavity_adapt

# Sources for executable
rayleigh_channel_SOURCES = rayleigh_channel.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rayleigh_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)



# Sources for executable
rayleigh_circular_driven_cavity_adapt_SOURCES = ray_circ_cavity_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rayleigh_circular_driven_cavity_adapt_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST+=flow.mcr
