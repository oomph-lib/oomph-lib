#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= \
driven_cavity \
direct_solver_test

# Do we have the TRILINOS library?
#=================================
if HAVE_TRILINOS
   check_PROGRAMS += TrilinosSolver_test
endif

#----------------------------------------------------------------------


# Sources for executable
driven_cavity_SOURCES = driven_cavity.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
driven_cavity_LDADD = \
                -L@libdir@ -lnavier_stokes  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Sources for executable
direct_solver_test_SOURCES = direct_solver_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
direct_solver_test_LDADD = \
                -L@libdir@ -lpoisson  \
                -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#----------------------------------------------------------------------


# Do we have the TRILINOS library?
#=================================
if HAVE_TRILINOS

# Sources the executable depends on:
TrilinosSolver_test_SOURCES = TrilinosSolver_test.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" library
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
TrilinosSolver_test_LDADD = -L@libdir@ -lpoisson -lgeneric \
                       $(EXTERNAL_LIBS) \
                       $(FLIBS)

endif



# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  

