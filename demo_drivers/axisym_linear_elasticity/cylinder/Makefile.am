#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=cylinder

# Sources for executable
cylinder_SOURCES = cylinder.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
cylinder_LDADD = -L@libdir@ -laxisym_linear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

