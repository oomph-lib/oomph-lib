#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=matrix_multiply_test 


# FIRST EXECUTABLE: SINGLE DRIVER FILE
#-------------------------------------


# Sources the executable depends on:
matrix_multiply_test_SOURCES = matrix_multiply_test.cc validate.sh

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources.
matrix_multiply_test_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
AM_CPPFLAGS += -I@includedir@

# Only run self test if we have python for fpdiff.py
#---------------------------------------------------
if HAVE_PYTHON
  # Test script
  TESTS =  validate.sh
else
  # Just run the executables
  TESTS = matrix_multiply_test
endif

