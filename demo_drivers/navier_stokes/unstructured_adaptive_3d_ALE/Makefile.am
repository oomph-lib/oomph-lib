#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=uns_adapt_3d

# Sources for executable
uns_adapt_3d_SOURCES = uns_adapt_3d.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
uns_adapt_3d_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                      -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

#EXTRA_DIST+=big_moving_block.mcr  check_navier_stokes_projection.mcr  \
#            moving_block.mcr
