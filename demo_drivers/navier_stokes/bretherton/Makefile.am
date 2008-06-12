#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=bretherton 

# Sources for executable
bretherton_SOURCES = bretherton.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
bretherton_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
