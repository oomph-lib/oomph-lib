#---------------------------------------------------------------------------
# Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Name of executable
check_PROGRAMS = space_time_oscillating_cylinder

#---------------------------------------------------------------------------

# Sources for executable
space_time_oscillating_cylinder_SOURCES = space_time_oscillating_cylinder.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
space_time_oscillating_cylinder_LDADD = -L@libdir@ \
-lspace_time_navier_stokes \
-lspace_time_block_preconditioner \
-lnavier_stokes \
-lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------
