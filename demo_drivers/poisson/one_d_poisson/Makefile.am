#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=one_d_poisson

# Sources for executable
one_d_poisson_SOURCES = one_d_poisson.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
