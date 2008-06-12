#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS=two_d_linear_wave_adapt

# Sources for executable
two_d_linear_wave_adapt_SOURCES = two_d_linear_wave_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_linear_wave_adapt_LDADD = -L@libdir@ -llinear_wave -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
