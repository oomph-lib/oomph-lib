#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
3d_ref_b_convect 3d_multimesh_ref_b_convect

#-------------------------------------------------------------------

# Sources for third executable
3d_ref_b_convect_SOURCES = 3d_ref_b_convect.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
3d_ref_b_convect_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------


# Sources for executable
3d_multimesh_ref_b_convect_SOURCES = 3d_multimesh_ref_b_convect.cc \
                                      my_boussinesq_elements.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
3d_multimesh_ref_b_convect_LDADD = -L@libdir@ -ladvection_diffusion \
  -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

