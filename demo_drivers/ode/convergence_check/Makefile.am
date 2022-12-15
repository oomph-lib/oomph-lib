#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=convergence_check

# Sources for executable
convergence_check_SOURCES = convergence_check.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
convergence_check_LDADD = -L@libdir@ -lode -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST+=zeros
