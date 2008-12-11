#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=pseudo_solid_collapsible_channel

#-------------------------------------------------------------------

# Sources for executable
pseudo_solid_collapsible_channel_SOURCES = pseudo_solid_collapsible_channel.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pseudo_solid_collapsible_channel_LDADD = -L@libdir@ -lnavier_stokes -lsolid  -lconstitutive -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------