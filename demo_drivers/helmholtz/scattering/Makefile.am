#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=scattering  adaptive_scattering 

#--------------------------------------

# Sources for executable
adaptive_scattering_SOURCES = scattering.cc helmholtz_bc_elements.h validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_scattering_LDADD = -L@libdir@ -lhelmholtz -lgeneric -loomph_crbond_bessel $(EXTERNAL_LIBS) $(FLIBS)

adaptive_scattering_CXXFLAGS = -DADAPTIVE



#--------------------------------------


# Sources for executable
scattering_SOURCES = scattering.cc helmholtz_bc_elements.h validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
scattering_LDADD = -L@libdir@ -lhelmholtz -lgeneric -loomph_crbond_bessel $(EXTERNAL_LIBS) $(FLIBS)


