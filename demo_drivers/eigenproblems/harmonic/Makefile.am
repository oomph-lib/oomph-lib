#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers


#=========================
# THESE ALL NEED ARPACK!
#=========================
if HAVE_ARPACK_SOURCES

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

endif

do_we_have_arpack:
if HAVE_ARPACK_SOURCES
	echo "have arpack"
else
	echo "dont have arpack"
endif
