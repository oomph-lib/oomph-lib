#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_three_d_solid  

#-------------------------------------------------------------------


# Sources for executable
unstructured_three_d_solid_SOURCES = unstructured_three_d_solid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_solid_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

EXTRA_DIST += three_d_unstructured_solid.mcr \
fsi_bifurcation_solid.1.ele  fsi_bifurcation_solid.1.face  \
fsi_bifurcation_solid.1.node  fsi_bifurcation_solid.poly


