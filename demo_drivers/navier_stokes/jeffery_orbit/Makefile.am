#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=jeffery_orbit

# Sources for executable
jeffery_orbit_SOURCES = jeffery_orbit.cc my_taylor_hood_elements.h


# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
jeffery_orbit_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                       -lrigid_body -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

#EXTRA_DIST+=big_moving_block.mcr  check_navier_stokes_projection.mcr  \
#            moving_block.mcr

