#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=\
oscillating_sphere \
adaptive_oscillating_sphere 


#--------------------------------------

# Sources for executable
oscillating_sphere_SOURCES = oscillating_sphere.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
oscillating_sphere_LDADD = -L@libdir@ -lpml_fourier_decomposed_helmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
adaptive_oscillating_sphere_SOURCES = oscillating_sphere.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_oscillating_sphere_LDADD = -L@libdir@ -lpml_fourier_decomposed_helmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_oscillating_sphere_CXXFLAGS = -DADAPTIVE

#--------------------------------------
#EXTRA_DIST += bessel.map exact.map

