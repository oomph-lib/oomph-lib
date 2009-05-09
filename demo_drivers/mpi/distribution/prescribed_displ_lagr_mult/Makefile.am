#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
prescribed_displ_lagr_mult

#-----------------------------------------------------------------------

# Sources for executable
prescribed_displ_lagr_mult_SOURCES = prescribed_displ_lagr_mult.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
prescribed_displ_lagr_mult_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------


EXTRA_DIST += presc_displ_lagr_mult_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  

