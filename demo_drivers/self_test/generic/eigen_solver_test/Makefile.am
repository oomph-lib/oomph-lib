#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= eigen_solver_test lapack_qz_find_eigenvalues_test 

#----------------------------------------------------------------------

# Sources for executable
eigen_solver_test_SOURCES = eigen_solver_test.cc
lapack_qz_find_eigenvalues_test_SOURCES = lapack_qz_find_eigenvalues_test.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
eigen_solver_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
lapack_qz_find_eigenvalues_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS += -DPARANOID
