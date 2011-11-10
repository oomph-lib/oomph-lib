#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=jeffery_orbit rot_free_surface

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


# Sources for executable
rot_free_surface_SOURCES = rot_free_surface.cc my_taylor_hood_elements.h


# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
rot_free_surface_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                       -lfluid_interface \
                       -lrigid_body -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)
