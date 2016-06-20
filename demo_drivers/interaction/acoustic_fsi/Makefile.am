#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS= \
acoustic_fsi \
unstructured_acoustic_fsi

#---------------------------------------------------------------------

# Sources for first executable
unstructured_acoustic_fsi_SOURCES = unstructured_acoustic_fsi.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_acoustic_fsi_LDADD =  -L@libdir@ -lhelmholtz \
                      -ltime_harmonic_linear_elasticity  \
                      -lnavier_stokes \
                      -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------------

# Sources for first executable
acoustic_fsi_SOURCES = acoustic_fsi.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
acoustic_fsi_LDADD =  -L@libdir@ -lhelmholtz \
                      -ltime_harmonic_linear_elasticity  \
                      -lnavier_stokes \
                      -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------

EXTRA_DIST += animate_displacement.mcr  animate_displacement_with_mesh.mcr \
 convert.bash 	run.bash animate_helmholtz.mcr

