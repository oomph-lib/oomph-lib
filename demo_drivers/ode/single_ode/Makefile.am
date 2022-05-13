#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=single_ode

# Sources for executable
single_ode_SOURCES = single_ode.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
single_ode_LDADD = -L@libdir@ -lode -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST+=zeros
