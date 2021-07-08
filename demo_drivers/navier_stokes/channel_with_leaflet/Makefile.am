#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=channel_with_leaflet 

# Sources for executable
channel_with_leaflet_SOURCES = channel_with_leaflet.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
channel_with_leaflet_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST += animate_channel_with_leaflet.mcr 