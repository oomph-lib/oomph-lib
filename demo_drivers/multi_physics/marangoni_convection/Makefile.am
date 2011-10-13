#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
marangoni_convection  marangoni_convection_box


#---------------------------------------------------------------------


# Sources for executable
marangoni_convection_SOURCES = marangoni_convection.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
marangoni_convection_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lfluid_interface -lmulti_physics -lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
marangoni_convection_box_SOURCES = marangoni_convection_box.cc 


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
marangoni_convection_box_LDADD = -L@libdir@ -ladvection_diffusion \
 -lnavier_stokes -lmulti_physics -lfluid_interface -lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)

