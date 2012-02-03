#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=hp_adaptive_driven_cavity

# Sources for executable
hp_adaptive_driven_cavity_SOURCES = hp_adaptive_driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
hp_adaptive_driven_cavity_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
