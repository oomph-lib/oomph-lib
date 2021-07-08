#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS=two_d_linear_wave two_d_linear_wave_flux

# Sources for first executable
two_d_linear_wave_SOURCES = two_d_linear_wave.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_linear_wave_LDADD = -L@libdir@ -llinear_wave -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for second executable
two_d_linear_wave_flux_SOURCES = two_d_linear_wave_flux.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_linear_wave_flux_LDADD = -L@libdir@ -llinear_wave -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


