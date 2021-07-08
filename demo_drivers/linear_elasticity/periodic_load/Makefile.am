#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=periodic_load refineable_periodic_load

# Sources for executable
periodic_load_SOURCES = periodic_load.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
periodic_load_LDADD = -L@libdir@ -llinear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
refineable_periodic_load_SOURCES = refineable_periodic_load.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
refineable_periodic_load_LDADD = -L@libdir@ -llinear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
