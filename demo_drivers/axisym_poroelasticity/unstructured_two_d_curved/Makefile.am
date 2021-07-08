#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
unstructured_two_d_curved \
adaptive_unstructured_two_d_curved 


#------------------------------------------------------


# Sources for executable
unstructured_two_d_curved_SOURCES = unstructured_two_d_curved.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_curved_LDADD = -L@libdir@ -laxisym_poroelasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------


# Sources for executable
adaptive_unstructured_two_d_curved_SOURCES = unstructured_two_d_curved.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_two_d_curved_LDADD = -L@libdir@ -laxisym_poroelasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_unstructured_two_d_curved_CXXFLAGS=-DADAPTIVE

#------------------------------------------------------

EXTRA_DIST += run_it.bash plot_it.bash \
animate.pvsm  validate.pvsm animate.py \
animate_validation_case.pvsm  \
validate_with_grouped_datasets.pvsm


