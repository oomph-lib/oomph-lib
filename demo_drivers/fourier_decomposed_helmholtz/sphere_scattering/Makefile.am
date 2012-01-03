#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=sphere_scattering 

#--------------------------------------

# Sources for executable
sphere_scattering_SOURCES = sphere_scattering.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
sphere_scattering_LDADD = -L@libdir@ -lfourier_decomposed_helmholtz -lgeneric -loomph_crbond_bessel $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

EXTRA_DIST += bessel.map