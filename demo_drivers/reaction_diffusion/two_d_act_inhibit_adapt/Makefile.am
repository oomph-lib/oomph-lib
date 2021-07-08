#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS= two_d_act_inhibit_adapt

# Sources for first executable
two_d_act_inhibit_adapt_SOURCES = two_d_act_inhibit_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_act_inhibit_adapt_LDADD = -L@libdir@ -ladvection_diffusion_reaction \
	                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
