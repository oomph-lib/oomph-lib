#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=helmholtz_point_source \
adaptive_helmholtz_point_source

#--------------------------------------

# Sources for executable
helmholtz_point_source_SOURCES = helmholtz_point_source.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
helmholtz_point_source_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS)  $(FLIBS)

#--------------------------------------

# Sources for executable
adaptive_helmholtz_point_source_SOURCES = helmholtz_point_source.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_helmholtz_point_source_LDADD = -L@libdir@ -lhelmholtz -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_helmholtz_point_source_CXXFLAGS = -DADAPTIVE

#--------------------------------------


