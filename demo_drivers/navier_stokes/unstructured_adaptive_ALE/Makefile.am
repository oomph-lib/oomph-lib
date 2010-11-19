#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_adaptive_ALE

# Sources for executable
unstructured_adaptive_ALE_SOURCES = unstructured_adaptive_ALE.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_adaptive_ALE_LDADD = -L@libdir@ \
                       -lsolid -lconstitutive -lnavier_stokes -lrigid_body \
                      -lgeneric \
                       $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST+=big_moving_block.mcr  check_navier_stokes_projection.mcr  \
            moving_block.mcr
