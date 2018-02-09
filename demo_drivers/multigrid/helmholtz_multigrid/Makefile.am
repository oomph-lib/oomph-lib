#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS= \
unstructured_two_d_helmholtz \
structured_cubic_point_source

#-----------------------------------------------------------------------------

# Sources for executable
unstructured_two_d_helmholtz_SOURCES = unstructured_two_d_helmholtz.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_helmholtz_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
structured_cubic_point_source_SOURCES = structured_cubic_point_source.cc

structured_cubic_point_source_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

#EXTRA_DIST += 
