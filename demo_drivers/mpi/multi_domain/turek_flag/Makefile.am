#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
turek_flag \
turek_flag_load_balance

#----------------------------------------------------------------------

# Sources for executable
turek_flag_SOURCES = turek_flag.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
turek_flag_LDADD = \
                -L@libdir@ -lsolid -lconstitutive -lnavier_stokes  \
                -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
turek_flag_load_balance_SOURCES = turek_flag_load_balance.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
turek_flag_load_balance_LDADD = \
                -L@libdir@ -lsolid -lconstitutive -lnavier_stokes  \
                -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

EXTRA_DIST += turek_flag_partition.dat load_balance.mcr

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS +=   -I@includedir@  

