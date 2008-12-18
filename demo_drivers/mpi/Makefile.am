#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= 

#=========================================
# THE WHOLE LOT IS ONLY RUN IF WE HAVE MPI
#=========================================
if WANT_MPI

# Do we have the TRILINOS library?
#=================================
if HAVE_TRILINOS
   check_PROGRAMS += TrilinosSolver_test
endif

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

endif


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES =   -I@includedir@  

