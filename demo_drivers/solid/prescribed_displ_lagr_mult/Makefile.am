#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=prescribed_displ_lagr_mult prescribed_displ_lagr_mult2  \
prescribed_displ_lagr_mult_precond

#-------------------------------------------------------------------

# Sources for executable
prescribed_displ_lagr_mult_SOURCES = prescribed_displ_lagr_mult.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
prescribed_displ_lagr_mult_LDADD = -L@libdir@ -lsolid -lconstitutive \
                          -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

# Sources for executable
prescribed_displ_lagr_mult2_SOURCES = prescribed_displ_lagr_mult2.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
prescribed_displ_lagr_mult2_LDADD = -L@libdir@ -lsolid -lconstitutive \
                       -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------------------------

# Sources for executable
prescribed_displ_lagr_mult_precond_SOURCES = prescribed_displ_lagr_mult_precond.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
prescribed_displ_lagr_mult_precond_LDADD = -L@libdir@ -lmulti_physics -lsolid -lconstitutive -lnavier_stokes -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

EXTRA_DIST+=lagr_distort.mcr run.sh
