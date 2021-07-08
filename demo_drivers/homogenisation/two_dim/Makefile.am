#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_dim two_dim_hex

# Sources for executable
two_dim_SOURCES = homo_lin_elasticity_elements.h \
homo_lin_elasticity_elements.cc  Thomo_lin_elasticity_elements.h \
Thomo_lin_elasticity_elements.cc two_dim.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_dim_LDADD = -L@libdir@ -llinear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
two_dim_hex_SOURCES = homo_lin_elasticity_elements.h \
homo_lin_elasticity_elements.cc  Thomo_lin_elasticity_elements.h \
Thomo_lin_elasticity_elements.cc two_dim_hex.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_dim_hex_LDADD = -L@libdir@ -llinear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

