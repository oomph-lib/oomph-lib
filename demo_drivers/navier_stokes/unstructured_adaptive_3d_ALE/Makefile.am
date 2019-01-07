#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=

if OOMPH_HAS_CGAL 

# Currently tet mesh adaptation only works with cgal based nearest
# sample point finder in solution transfer

check_PROGRAMS+=\
uns_adapt_3d # uns_adapt_3d_fs

endif

# Sources for executable
uns_adapt_3d_SOURCES = uns_adapt_3d.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
uns_adapt_3d_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes \
                      -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
#uns_adapt_3d_fs_SOURCES = uns_adapt_3d_fs.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
#uns_adapt_3d_fs_LDADD = -L@libdir@ \
#                       -lsolid -lconstitutive -lnavier_stokes \
#                      -lfluid_interface -lgeneric \
#                       $(EXTERNAL_LIBS) $(FLIBS)

#EXTRA_DIST+=big_moving_block.mcr  check_navier_stokes_projection.mcr  \
#            moving_block.mcr
