#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_two_d_helmholtz \
unstructured_two_d_helmholtz_adapt \
unstructured_two_d_helmholtz_scattering

#--------------------------------------

# Sources for executable
unstructured_two_d_helmholtz_scattering_SOURCES =  unstructured_two_d_helmholtz_scattering.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_helmholtz_scattering_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
unstructured_two_d_helmholtz_SOURCES =  unstructured_two_d_helmholtz.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_helmholtz_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
unstructured_two_d_helmholtz_adapt_SOURCES =  unstructured_two_d_helmholtz.cc


# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_two_d_helmholtz_adapt_LDADD = -L@libdir@ -lpml_helmholtz -lmeshes -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

unstructured_two_d_helmholtz_adapt_CXXFLAGS = -DADAPTIVE


#--------------------------------------


#EXTRA_DIST += annular_meshes.h
