#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=eighth_sphere_poisson

# Sources for executable
eighth_sphere_poisson_SOURCES = eighth_sphere_poisson.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
eighth_sphere_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric \
                              $(EXTERNAL_LIBS) $(FLIBS)
