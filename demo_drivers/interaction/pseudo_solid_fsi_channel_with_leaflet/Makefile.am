#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
fsi_channel_with_leaflet_precond

#---------------------------------------------------------------------

# Sources for first executable
fsi_channel_with_leaflet_precond_SOURCES = fsi_channel_with_leaflet_precond.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_channel_with_leaflet_precond_LDADD = \
 -L@libdir@ -lmulti_physics -lbeam -lnavier_stokes -lconstitutive -lsolid -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------

EXTRA_DIST += animate_channel_with_leaflet.mcr run.sh
