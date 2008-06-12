#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS= two_d_adv_diff_adapt two_d_adv_diff_adapt2

# Sources for first executable
two_d_adv_diff_adapt_SOURCES = two_d_adv_diff_adapt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_adv_diff_adapt_LDADD = -L@libdir@ -ladvection_diffusion -lgeneric \
                             $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
two_d_adv_diff_adapt2_SOURCES = two_d_adv_diff_adapt2.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_adv_diff_adapt2_LDADD = -L@libdir@ -ladvection_diffusion -lgeneric \
                              $(EXTERNAL_LIBS) $(FLIBS)
