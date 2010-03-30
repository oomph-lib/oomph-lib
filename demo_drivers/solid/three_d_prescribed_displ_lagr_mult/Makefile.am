#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=three_d_prescribed_displ_lagr_mult  \
three_d_prescribed_displ_lagr_mult_non_adapt 

#-------------------------------------------------------------------

# Sources for executable
three_d_prescribed_displ_lagr_mult_SOURCES = three_d_prescribed_displ_lagr_mult.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_prescribed_displ_lagr_mult_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------


# Sources for executable
three_d_prescribed_displ_lagr_mult_non_adapt_SOURCES = three_d_prescribed_displ_lagr_mult.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_prescribed_displ_lagr_mult_non_adapt_LDADD = -L@libdir@ -lsolid -lconstitutive -lgeneric \
                    $(EXTERNAL_LIBS) $(FLIBS)


# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_prescribed_displ_lagr_mult_non_adapt_CXXFLAGS = -DNOREFINE

#-------------------------------------------------------------------

EXTRA_DIST+=block_prescr.mcr
