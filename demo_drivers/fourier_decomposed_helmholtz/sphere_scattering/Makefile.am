#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=sphere_scattering  unstructured_sphere_scattering  \
unstructured_adaptive_sphere_scattering 

#--------------------------------------

# Sources for executable
sphere_scattering_SOURCES = sphere_scattering.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
sphere_scattering_LDADD = -L@libdir@ -lfourier_decomposed_helmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
unstructured_sphere_scattering_SOURCES = unstructured_sphere_scattering.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_sphere_scattering_LDADD = -L@libdir@ -lfourier_decomposed_helmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#--------------------------------------

# Sources for executable
unstructured_adaptive_sphere_scattering_SOURCES = unstructured_sphere_scattering.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_adaptive_sphere_scattering_LDADD = -L@libdir@ -lfourier_decomposed_helmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

unstructured_adaptive_sphere_scattering_CXXFLAGS = -DADAPTIVE

#--------------------------------------

EXTRA_DIST += bessel.map exact.map

