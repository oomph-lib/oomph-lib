#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=two_d_tilted_square


#---------------------------------------------------------------

# Sources for first executable
#two_d_fp_SOURCES = two_d_fp.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
#                     -lnavier_stokes -lgeneric \
#                     $(EXTERNAL_LIBS)  $(FLIBS)

#---------------------------------------------------------------

# Sources for first executable
two_d_tilted_square_SOURCES = two_d_tilted_square.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_tilted_square_LDADD = -L@libdir@ \
                     -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS)  $(FLIBS)

#---------------------------------------------------------------


EXTRA_DIST+=run.sh


#EXTRA_DIST+=analyse_fp_2d.mcr  analyse_fp_3d.mcr  animate_3d.mcr


