#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


# Name of executable
check_PROGRAMS=harmonic complex_harmonic

# Sources for executable
harmonic_SOURCES = harmonic.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
harmonic_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
complex_harmonic_SOURCES = complex_harmonic.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
complex_harmonic_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
