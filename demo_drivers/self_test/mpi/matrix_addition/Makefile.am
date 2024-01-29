#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
matrix_addition

#----------------------------------------------------------------------

# Sources for executable
matrix_addition_SOURCES = matrix_addition.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
matrix_addition_LDADD = \
                -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

