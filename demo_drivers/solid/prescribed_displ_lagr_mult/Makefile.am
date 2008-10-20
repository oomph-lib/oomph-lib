#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=prescribed_displ_lagr_mult

#-------------------------------------------------------------------

# Sources for executable
prescribed_displ_lagr_mult_SOURCES = prescribed_displ_lagr_mult.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
prescribed_displ_lagr_mult_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

EXTRA_DIST+=lagr_distort.mcr
