#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=adaptive_fs adaptive_two_fluid

# Sources for executable
adaptive_fs_SOURCES = adaptive_fs.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_fs_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                       -lfluid_interface -lrigid_body -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
adaptive_two_fluid_SOURCES = adaptive_two_fluid.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_two_fluid_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                       -lfluid_interface -lrigid_body -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)


#EXTRA_DIST+=big_moving_block.mcr  check_navier_stokes_projection.mcr  \
#            moving_block.mcr
