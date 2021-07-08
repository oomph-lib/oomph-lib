#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= unstructured_three_d_fsi


#---------------------------------------------------------------------

# Sources for first executable
unstructured_three_d_fsi_SOURCES = unstructured_three_d_fsi.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_three_d_fsi_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#unstructured_three_d_fsi_CXXFLAGS= -g -DPARANOID -DRANGE_CHECKING

#---------------------------------------------------------------------


EXTRA_DIST += \
fsi_bifurcation_fluid.1.ele   fsi_bifurcation_fluid.1.node  \
fsi_bifurcation_solid.1.face  fsi_bifurcation_solid.poly \
fsi_bifurcation_fluid.1.face  fsi_bifurcation_fluid.poly  \
fsi_bifurcation_solid.1.ele  fsi_bifurcation_solid.1.node 




