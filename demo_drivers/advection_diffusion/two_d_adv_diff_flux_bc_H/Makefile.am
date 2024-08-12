#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= two_d_adv_diff_flux_bc

# Sources for executable
two_d_adv_diff_flux_bc_SOURCES = two_d_adv_diff_flux_bc.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_adv_diff_flux_bc_LDADD = -L@libdir@ -ladvection_diffusion -lgeneric \
                               $(EXTERNAL_LIBS) $(FLIBS)
