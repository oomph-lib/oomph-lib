#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=jeffery_hamel

# Sources for executable
jeffery_hamel_SOURCES = jeffery_hamel.cc refineable_r_mesh.h \
refineable_streamfunction_elements.h streamfunction_elements.h \
streamfunction_include.h jh_includes.h polar_streamfunction_traction_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
jeffery_hamel_LDADD = -L@libdir@ -lpolar_navier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
