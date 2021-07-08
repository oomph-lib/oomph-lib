#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
three_d_cantilever_adapt

#----------------------------------------------------------------------

# Sources for executable
three_d_cantilever_adapt_SOURCES = three_d_cantilever.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_cantilever_adapt_LDADD = \
                -L@libdir@ -lsolid -lconstitutive  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
three_d_cantilever_adapt_CXXFLAGS = -DREFINE

#----------------------------------------------------------------------

EXTRA_DIST += three_d_cantilever_0_partition.dat \
              three_d_cantilever_1_partition.dat \
              three_d_cantilever_3_partition.dat \
              three_d_cantilever_4_partition.dat \
              three_d_cantilever_5_partition.dat \
              three_d_cantilever_6_partition.dat \
              three_d_cantilever_8_partition.dat \
              three_d_cantilever_9_partition.dat \
              comp.mcr
# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

