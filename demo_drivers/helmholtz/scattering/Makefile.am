#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=scattering  adaptive_scattering unstructured_scattering \
adaptive_unstructured_scattering

#--------------------------------------

# Sources for executable
adaptive_scattering_SOURCES = scattering.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_scattering_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_scattering_CXXFLAGS = -DADAPTIVE

#--------------------------------------

# Sources for executable
scattering_SOURCES = scattering.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
scattering_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
unstructured_scattering_SOURCES = unstructured_scattering.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_scattering_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
adaptive_unstructured_scattering_SOURCES = unstructured_scattering.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_scattering_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_unstructured_scattering_CXXFLAGS = -DADAPTIVE


#--------------------------------------


