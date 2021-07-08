#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS= two_d_variable_diff_adapt

# Sources for first executable
two_d_variable_diff_adapt_SOURCES = two_d_variable_diff_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_variable_diff_adapt_LDADD = -L@libdir@ -ladvection_diffusion -lgeneric \
                             $(EXTERNAL_LIBS) $(FLIBS)

