#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= unstructured_two_d_fsi


#---------------------------------------------------------------------

# Sources for first executable
unstructured_two_d_fsi_SOURCES = unstructured_two_d_fsi.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_fsi_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------


EXTRA_DIST += fluid.fig  fluid.fig.1.ele  fluid.fig.1.node  \
fluid.fig.1.poly	fluid.fig.poly	solid.fig  solid.fig.1.ele  \
solid.fig.1.node  solid.fig.1.poly  solid.fig.poly
