#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_d_biharmonic 

# Sources for executable
two_d_biharmonic_SOURCES = two_d_biharmonic.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_biharmonic_LDADD = -L@libdir@ -lbiharmonic -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
