#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=one_d_spectral two_d_spectral_adapt three_d_spectral_adapt

# Sources for executable
one_d_spectral_SOURCES = one_d_spectral.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_spectral_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
two_d_spectral_adapt_SOURCES = two_d_spectral_adapt.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_spectral_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
three_d_spectral_adapt_SOURCES = three_d_spectral_adapt.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_spectral_adapt_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
