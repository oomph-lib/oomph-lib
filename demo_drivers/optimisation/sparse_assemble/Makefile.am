#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=sparse_assemble_test

# Sources for executable
sparse_assemble_test_SOURCES = sparse_assemble_test.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
sparse_assemble_test_LDADD = -L@libdir@ -lnavier_stokes \
	        	     -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

EXTRA_DIST+=compare.bash


