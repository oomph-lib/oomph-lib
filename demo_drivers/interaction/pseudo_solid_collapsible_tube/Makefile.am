#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
pseudo_solid_collapsible_tube \
pseudo_solid_collapsible_tube_work


#---------------------------------------------------------------------

# Sources for executable
pseudo_solid_collapsible_tube_SOURCES = pseudo_solid_collapsible_tube.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pseudo_solid_collapsible_tube_LDADD = -L@libdir@ -lmulti_physics -lnavier_stokes -lsolid  -lconstitutive -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

# Sources for executable
pseudo_solid_collapsible_tube_work_SOURCES = pseudo_solid_collapsible_tube_work.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pseudo_solid_collapsible_tube_work_LDADD = -L@libdir@ -lmulti_physics -lnavier_stokes -lsolid  -lconstitutive -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------

EXTRA_DIST += pseudo_solid_collapsible_tube.mcr
