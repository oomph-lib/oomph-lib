#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=rotating_ends 

# Sources for executable
rotating_ends_SOURCES = rotating_ends.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
rotating_ends_LDADD = -L@libdir@ -laxisym_navier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)
