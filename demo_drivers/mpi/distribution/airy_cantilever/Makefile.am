#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
airy_cantilever2_adapt

#----------------------------------------------------------------------

# Sources for executable
airy_cantilever2_adapt_SOURCES = airy_cantilever2.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever2_adapt_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
airy_cantilever2_adapt_CXXFLAGS = -DREFINE

#-----------------------------------------------------------------------

EXTRA_DIST += airy_cantilever_0_partition.dat \
              airy_cantilever_1_partition.dat \
              airy_cantilever_2_partition.dat \
              airy_cantilever_3_partition.dat \
              airy_cantilever_4_partition.dat \
              airy_cantilever_5_partition.dat \
              airy_cantilever_6_partition.dat \
              airy_cantilever_7_partition.dat \
              airy_cantilever_8_partition.dat \
              airy_cantilever_9_partition.dat

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

