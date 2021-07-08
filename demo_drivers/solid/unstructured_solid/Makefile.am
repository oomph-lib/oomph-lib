#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_three_d_solid  unstructured_two_d_solid


#-------------------------------------------------------------------

# Sources for executable
unstructured_two_d_solid_SOURCES = unstructured_two_d_solid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_solid_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------


# Sources for executable
unstructured_three_d_solid_SOURCES = unstructured_three_d_solid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_solid_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

EXTRA_DIST += solid.fig  solid.fig.1.ele  solid.fig.1.node  \
solid.fig.1.poly solid.fig.poly cube_hole.1.ele  \
cube_hole.1.face  cube_hole.1.node \
three_d_unstructured_solid.mcr  two_d_unstructured_solid.mcr

