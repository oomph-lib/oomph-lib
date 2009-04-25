#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_three_d_fluid

# Sources for executable
unstructured_three_d_fluid_SOURCES = unstructured_three_d_fluid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_fluid_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST+= fsi_bifurcation_fluid.poly  \
fsi_bifurcation_fluid.1.ele  \
fsi_bifurcation_fluid.1.node \
fsi_bifurcation_fluid.1.face \
fsi_bifurcation_fluid_with_extended_tubes.poly   \
fsi_bifurcation_fluid_with_extended_tubes.1.node \
fsi_bifurcation_fluid_with_extended_tubes.1.ele  \
fsi_bifurcation_fluid_with_extended_tubes.1.face \
three_d_unstructured_fluid.mcr
