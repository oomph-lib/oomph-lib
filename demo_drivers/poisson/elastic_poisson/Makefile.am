#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=elastic_poisson 

# Sources for executable
elastic_poisson_SOURCES = elastic_poisson.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
elastic_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
