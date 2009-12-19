#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


# Programs with self tests
check_PROGRAMS= steady_axisym_advection_diffusion

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources


#Sources
steady_axisym_advection_diffusion_SOURCES =\
steady_axisym_advection_diffusion.cc

# Required libraries: 
steady_axisym_advection_diffusion_LDADD = -L@libdir@ \
                                          -lsteady_axisym_advection_diffusion \
                                          -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

