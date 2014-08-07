#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_layer_interface_nonaxisym_perturbations

# Sources for executable
two_layer_interface_nonaxisym_perturbations_SOURCES = \
two_layer_interface_nonaxisym_perturbations.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_layer_interface_nonaxisym_perturbations_LDADD = -L@libdir@ \
 -laxisym_navier_stokes -lnavier_stokes -lfluid_interface \
 -lgeneric -loomph_crbond_bessel $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
INCLUDES = -I@includedir@  


