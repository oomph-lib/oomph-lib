#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS = \
unit_cube_poisson \
unit_cube_poisson_validate \
two_d_poisson_tanh_flux_bc \
two_d_poisson_tanh_flux_bc_validate

#---------------------------------------------------------------------------

# Sources for executable
unit_cube_poisson_SOURCES = unit_cube_poisson.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unit_cube_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for executable
unit_cube_poisson_validate_SOURCES = unit_cube_poisson_validate.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unit_cube_poisson_validate_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Local sources that the MG driver code depends on:
two_d_poisson_tanh_flux_bc_SOURCES = two_d_poisson_tanh_flux_bc.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson_tanh_flux_bc_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Local sources that the MG driver code depends on:
two_d_poisson_tanh_flux_bc_validate_SOURCES = two_d_poisson_tanh_flux_bc_validate.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson_tanh_flux_bc_validate_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------
