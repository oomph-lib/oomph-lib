#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executables
if OOMPH_HAS_HYPRE
if OOMPH_HAS_TRILINOS

check_PROGRAMS= \
pseudo_solid_collapsible_tube 

endif
endif

#---------------------------------------------------------------------

# Sources for executable
pseudo_solid_collapsible_tube_SOURCES = pseudo_solid_collapsible_tube.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
pseudo_solid_collapsible_tube_LDADD = -L@libdir@ -lmulti_physics -lnavier_stokes -lsolid  \
-lconstitutive -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------

EXTRA_DIST += pseudo_solid_collapsible_tube.mcr make_dummy_solid_solns.bash



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

