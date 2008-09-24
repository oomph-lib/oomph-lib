#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=spherical_couette

# Sources for executable
spherical_couette_SOURCES = spherical_couette.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spherical_couette_LDADD = -L@libdir@ -lspherical_navier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)
