#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=driven_cavity 

# Sources for executable
driven_cavity_SOURCES = driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
