#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=vmtk_fluid

# Sources for executable
vmtk_fluid_SOURCES = vmtk_fluid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
vmtk_fluid_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST+= fluid_iliac.1.ele  fluid_iliac.1.face  fluid_iliac.1.node  fluid_iliac.poly iliac.xda


