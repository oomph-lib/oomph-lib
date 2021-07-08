#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=time_harmonic_elastic_annulus \
adaptive_time_harmonic_elastic_annulus \
unstructured_time_harmonic_elastic_annulus \
adaptive_unstructured_time_harmonic_elastic_annulus 

#--------------------------------------

# Sources for executable
time_harmonic_elastic_annulus_SOURCES = time_harmonic_elastic_annulus.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
time_harmonic_elastic_annulus_LDADD = -L@libdir@ -ltime_harmonic_linear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------

# Sources for executable
adaptive_time_harmonic_elastic_annulus_SOURCES = time_harmonic_elastic_annulus.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_time_harmonic_elastic_annulus_LDADD = -L@libdir@ -ltime_harmonic_linear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

adaptive_time_harmonic_elastic_annulus_CXXFLAGS = -DADAPTIVE

#--------------------------------------


# Sources for executable
unstructured_time_harmonic_elastic_annulus_SOURCES = unstructured_time_harmonic_elastic_annulus.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_time_harmonic_elastic_annulus_LDADD = -L@libdir@ -ltime_harmonic_linear_elasticity -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


#--------------------------------------


# Sources for executable
adaptive_unstructured_time_harmonic_elastic_annulus_SOURCES = unstructured_time_harmonic_elastic_annulus.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
adaptive_unstructured_time_harmonic_elastic_annulus_LDADD = -L@libdir@ -ltime_harmonic_linear_elasticity -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

adaptive_unstructured_time_harmonic_elastic_annulus_CXXFLAGS = -DADAPTIVE

#--------------------------------------


EXTRA_DIST += animate_displacement.mcr elast.map maple.dat \
              animate_displacement_with_mesh.mcr convert.bash
