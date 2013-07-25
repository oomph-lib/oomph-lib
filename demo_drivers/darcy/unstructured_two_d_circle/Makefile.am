#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_two_d_circle adaptive_unstructured_two_d_circle

# Sources for executable
unstructured_two_d_circle_SOURCES = unstructured_two_d_circle.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_circle_LDADD = -L@libdir@ -ldarcy -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
adaptive_unstructured_two_d_circle_SOURCES = unstructured_two_d_circle.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_two_d_circle_LDADD = -L@libdir@ -ldarcy -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_unstructured_two_d_circle_CXXFLAGS=-DADAPTIVE


EXTRA_DIST += continuity_of_normal_flux_checker.pvsm p_shape.pvsm \
q_shape.pvsm plot_shape_fcts.bash flow.pvsm


