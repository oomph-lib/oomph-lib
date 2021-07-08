#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=rayleigh_channel

# Sources for executable
rayleigh_channel_SOURCES = rayleigh_channel.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rayleigh_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

