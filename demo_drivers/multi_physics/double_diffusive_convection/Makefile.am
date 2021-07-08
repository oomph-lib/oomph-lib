#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= dd_convection multimesh_dd_convection \
 multimesh_ref_dd_convection

#---------------------------------------------------------------------

# Sources for first executable
dd_convection_SOURCES = dd_convection.cc db_navier_st_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
dd_convection_LDADD = -L@libdir@ -ladvection_diffusion_reaction \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for first executable
multimesh_dd_convection_SOURCES = multimesh_dd_convection.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multimesh_dd_convection_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for first executable
multimesh_ref_dd_convection_SOURCES = multimesh_ref_dd_convection.cc \
      db_nst_external_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
multimesh_ref_dd_convection_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
