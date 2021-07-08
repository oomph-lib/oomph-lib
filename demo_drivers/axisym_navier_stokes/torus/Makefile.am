#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=torus

# Sources for executable
torus_SOURCES = torus.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
torus_LDADD = -L@libdir@ -lnavier_stokes -laxisym_navier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)
