#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS= two_d_adv_diff_SUPG

# Sources for first executable
two_d_adv_diff_SUPG_SOURCES = two_d_adv_diff_SUPG.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_adv_diff_SUPG_LDADD = -L@libdir@ -ladvection_diffusion -lgeneric \
                             $(EXTERNAL_LIBS) $(FLIBS)
