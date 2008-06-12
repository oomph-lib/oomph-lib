#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=one_d_poisson_generic_only

# Sources for executable
one_d_poisson_generic_only_SOURCES = one_d_poisson_generic_only.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_poisson_generic_only_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
