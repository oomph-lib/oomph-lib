#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable:
check_PROGRAMS= biharmonic_problem_test 

# Sources for executable
biharmonic_problem_test_SOURCES = \
biharmonic_problem_test.cc \
hermite_element_quad_mesh.h \
one_d_hermite_element_mesh.h  \
topologically_rectangular_domain.h


# Required libraries: 
biharmonic_problem_test_LDADD = -L@libdir@ -lbiharmonic -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
