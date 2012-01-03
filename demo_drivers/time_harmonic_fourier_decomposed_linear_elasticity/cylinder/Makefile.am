#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
cylinder pressure_loaded_cylinder


#--------------------------------------

# Sources for executable
cylinder_SOURCES = cylinder.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
cylinder_LDADD = -L@libdir@ -ltime_harmonic_fourier_decomposed_linear_elasticity -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
pressure_loaded_cylinder_SOURCES = pressure_loaded_cylinder.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pressure_loaded_cylinder_LDADD = -L@libdir@ -ltime_harmonic_fourier_decomposed_linear_elasticity -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# EXTRA_DIST += animate_displacement.mcr
