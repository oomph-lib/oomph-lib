#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=fish_poisson fish_poisson_no_adapt

# Sources for first executable
fish_poisson_SOURCES = fish_poisson.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
fish_poisson_no_adapt_SOURCES = fish_poisson_no_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fish_poisson_no_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric \
                              $(EXTERNAL_LIBS) $(FLIBS)
