#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# DO NOT NEED TO CHECK FOR MPI BECAUSE IF WE DO NOT HAVE MPI WE DO NOT
# DESCEND INTO THIS DIRECTORY

# Name of executable
check_PROGRAMS= eigen_solver_test lapack_qz_find_eigenvalues_test anasazi_mwe harmonic

#----------------------------------------------------------------------

# Sources for executable
eigen_solver_test_SOURCES = eigen_solver_test.cc
lapack_qz_find_eigenvalues_test_SOURCES = lapack_qz_find_eigenvalues_test.cc
anasazi_mwe_SOURCES = anasazi_mwe.cc
harmonic_SOURCES = harmonic.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
eigen_solver_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
lapack_qz_find_eigenvalues_test_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
anasazi_mwe_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)
harmonic_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) \
								  $(FLIBS)

#----------------------------------------------------------------------

# Include path for library headers: All library headers live in
# the include directory which we specify with -I
AM_CPPFLAGS += -DPARANOID
