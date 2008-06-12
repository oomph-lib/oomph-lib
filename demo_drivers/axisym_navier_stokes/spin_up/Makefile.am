#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=spin_up 

# Sources for executable
spin_up_SOURCES = spin_up.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
spin_up_LDADD = -L@libdir@ -laxisym_navier_stokes -lnavier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)
