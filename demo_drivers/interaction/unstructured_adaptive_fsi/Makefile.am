#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= unstructured_adaptive_2d_fsi


#---------------------------------------------------------------------

# Sources for first executable
unstructured_adaptive_2d_fsi_SOURCES = unstructured_adaptive_2d_fsi.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_adaptive_2d_fsi_LDADD = \
 -L@libdir@ -lsolid -lconstitutive -lnavier_stokes -lmeshes -lgeneric  \
                     $(EXTERNAL_LIBS) $(FLIBS)

