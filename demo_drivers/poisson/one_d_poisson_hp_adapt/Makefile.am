#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= one_d_poisson_hp_adapt

# Sources for executable
one_d_poisson_hp_adapt_SOURCES = one_d_poisson_hp_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_poisson_hp_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
