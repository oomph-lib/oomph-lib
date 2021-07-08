#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables 
check_PROGRAMS=two_d_poisson two_d_poisson2 two_d_poisson_compare_solvers

# Sources for first executable
two_d_poisson_SOURCES = two_d_poisson.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for second executable
two_d_poisson2_SOURCES = two_d_poisson2.cc two_d_poisson2_mesh.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson2_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for third executable 
two_d_poisson_compare_solvers_SOURCES = two_d_poisson_compare_solvers.cc 

# Do we have the hsl sources? 
if HAVE_HSL_SOURCES
  AM_CPPFLAGS+="-DHAVE_HSL_SOURCES"
endif

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_poisson_compare_solvers_LDADD = -L@libdir@ -lpoisson -lgeneric \
                                      $(EXTERNAL_LIBS) $(FLIBS)

