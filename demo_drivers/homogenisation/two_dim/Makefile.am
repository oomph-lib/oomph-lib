#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_dim

# Sources for executable
two_dim_SOURCES = homo_lin_elasticity_elements.h \
homo_lin_elasticity_elements.cc  Thomo_lin_elasticity_elements.h \
Thomo_lin_elasticity_elements.cc homo_lin_elasticity_traction_elemetns.h \
two_dim.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_dim_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

