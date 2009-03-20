#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_two_d_fluid

# Sources for executable
unstructured_two_d_fluid_SOURCES = unstructured_two_d_fluid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_fluid_LDADD = -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST+=fluid.fig  fluid.fig.1.ele  fluid.fig.1.node  \
            fluid.fig.1.poly	fluid.fig.poly
