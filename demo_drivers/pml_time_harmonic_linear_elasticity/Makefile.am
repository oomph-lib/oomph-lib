#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=time_harmonic_elasticity_driver_source \
time_harmonic_elasticity_driver

#--------------------------------------

# Sources for executable
time_harmonic_elasticity_driver_source_SOURCES = time_harmonic_elasticity_driver_source.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
time_harmonic_elasticity_driver_source_LDADD = -L@libdir@ -lpml_time_harmonic_linear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
time_harmonic_elasticity_driver_SOURCES = time_harmonic_elasticity_driver.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
time_harmonic_elasticity_driver_LDADD = -L@libdir@ -lpml_time_harmonic_linear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# EXTRA_DIST += animate_displacement.mcr elast.map maple.dat 
