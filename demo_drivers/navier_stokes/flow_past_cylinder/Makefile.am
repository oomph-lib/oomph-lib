#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=flow_past_cylinder

# Sources for executable
flow_past_cylinder_SOURCES = flow_past_cylinder.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
flow_past_cylinder_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

EXTRADIST = flow.lay