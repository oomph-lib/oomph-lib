#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


# Programs with self tests
check_PROGRAMS= thermal_fibre

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources


#Sources
thermal_fibre_SOURCES =\
non_isothermal_axisym_nst_elements.h thermal_fibre.cc

# Required libraries: 
thermal_fibre_LDADD = -L@libdir@ \
		      -laxisym_navier_stokes \
		      -lfluid_interface \
	              -lsteady_axisym_advection_diffusion \
	              -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

